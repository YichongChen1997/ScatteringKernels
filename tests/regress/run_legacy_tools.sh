#!/usr/bin/env bash
# Runs every legacy tool on one stored reference set, the same way the
# reference outputs were made (see tests/reference/R3_channel_mini/make.sh).
#
# Usage: run_legacy_tools.sh REFDIR BINDIR WORKDIR
#   REFDIR   a reference "mini" directory, e.g. tests/reference/R3_channel_mini/mini,
#            with DUMP_CONTENT.sha256, the gzipped dumps and Specification.dat
#   BINDIR   directory with the compiled tools (make puts them in build/bin)
#   WORKDIR  where to run. Only the files this script writes are replaced.
# Environment:
#   SKIP_TOOLS    comma-separated tool names not to run (default: none)
#   SKIP_MISSING_INPUT=1  do not run pp_meas_wallTemp when the set has no
#                 wall dump. Without its input the tool works on values
#                 that were never set, which change with the compiler
#                 (default 0: run it anyway, as the reference was made)
#   TOOL_TIMEOUT  seconds before a tool is killed (default 600)
#
# What it writes, laid out like the stored outputs/ so that the two can be
# compared file by file (tests/regress/compare_outputs.py):
#   WORKDIR/dump_meas_gas.lammpstrj    gas dump, unpacked and checked against
#   WORKDIR/dump_meas_wall.lammpstrj   DUMP_CONTENT.sha256 (wall dump if stored)
#   WORKDIR/Specification.dat          copied from REFDIR
#   WORKDIR/tools/<tool>/              the tool's own output files, _stdout.txt,
#                                      _stderr.txt (only if not empty), and links
#                                      to the three input files above
#   WORKDIR/tool_status.tsv            tool, exit_code, wall_s, note (a tool
#                                      that was not run has exit_code NA)
#
# Each tool runs in its own directory, with stdin from /dev/null. Needs bash
# (3.2 or later), perl, gzip, and sha256sum or shasum.
set -eu
export LC_ALL=C

usage() { sed -n '2,30p' "$0" >&2; exit 2; }
[ $# -eq 3 ] || usage

REF=$(cd "$1" && pwd)
BIN=$(cd "$2" && pwd)
mkdir -p "$3"
WORK=$(cd "$3" && pwd)
SKIP_TOOLS=${SKIP_TOOLS:-}
SKIP_MISSING_INPUT=${SKIP_MISSING_INPUT:-0}
TOOL_TIMEOUT=${TOOL_TIMEOUT:-600}

GAS=dump_meas_gas.lammpstrj
WALL=dump_meas_wall.lammpstrj
SPEC=Specification.dat

die() { echo "run_legacy_tools.sh: $*" >&2; exit 1; }
now() { perl -MTime::HiRes=time -e 'printf "%.3f\n", time'; }
sha_of() {
    if command -v sha256sum > /dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'
    else shasum -a 256 "$1" | awk '{print $1}'; fi
}
is_skipped() {
    case ",$SKIP_TOOLS," in *",$1,"*) return 0 ;; *) return 1 ;; esac
}

[ -f "$REF/DUMP_CONTENT.sha256" ] || die "no DUMP_CONTENT.sha256 in $REF"
[ -f "$REF/$SPEC" ] || die "no $SPEC in $REF"

# Clear only what this script writes.
rm -rf "$WORK/tools"
rm -f "$WORK/$GAS" "$WORK/$WALL" "$WORK/$SPEC" "$WORK/tool_status.tsv"
mkdir -p "$WORK/tools"

# Unpack the dumps under the names the tools open, and check their content.
while read -r sum name; do
    [ -n "$sum" ] || continue
    case "$name" in
        *wall*) dest=$WALL ;;
        *gas*)  dest=$GAS ;;
        *) die "cannot tell whether $name is the gas or the wall dump" ;;
    esac
    if [ -f "$REF/$name.gz" ]; then
        gzip -dc "$REF/$name.gz" > "$WORK/$dest"
    elif [ -f "$REF/$name" ]; then
        cp "$REF/$name" "$WORK/$dest"
    else
        die "neither $name.gz nor $name is in $REF"
    fi
    got=$(sha_of "$WORK/$dest")
    [ "$got" = "$sum" ] || die "$dest: sha256 $got, expected $sum (from $name)"
done < "$REF/DUMP_CONTENT.sha256"
cp "$REF/$SPEC" "$WORK/$SPEC"

# The tools to run: those listed in the reference tool_status.tsv, otherwise
# every pp_meas_* program in BINDIR.
if [ -f "$REF/tool_status.tsv" ]; then
    TOOLS=$(awk -F'\t' 'NR > 1 && $1 != "" {print $1}' "$REF/tool_status.tsv")
else
    TOOLS=$(cd "$BIN" && ls pp_meas_* 2> /dev/null || true)
fi
[ -n "$TOOLS" ] || die "no tools to run"

status=$WORK/tool_status.tsv
printf 'tool\texit_code\twall_s\tnote\n' > "$status"
for name in $TOOLS; do
    if is_skipped "$name"; then
        printf '%s\tNA\tNA\tskipped\n' "$name" >> "$status"
        echo "  $name skipped"
        continue
    fi
    if [ "$SKIP_MISSING_INPUT" = 1 ] && [ "$name" = pp_meas_wallTemp ] && [ ! -e "$WORK/$WALL" ]; then
        printf '%s\tNA\tNA\tskipped: no wall dump\n' "$name" >> "$status"
        echo "  $name skipped (no wall dump)"
        continue
    fi
    if [ ! -x "$BIN/$name" ]; then
        printf '%s\tNA\tNA\tnot built\n' "$name" >> "$status"
        echo "  $name not built"
        continue
    fi
    wd=$WORK/tools/$name
    mkdir -p "$wd"
    for f in "$GAS" "$WALL" "$SPEC"; do
        if [ -e "$WORK/$f" ]; then ln -s "../../$f" "$wd/$f"; fi
    done
    rc=0
    t0=$(now)
    # perl sets an alarm that survives exec, so a stuck tool is killed with
    # SIGALRM after TOOL_TIMEOUT seconds (exit code 142). Waiting for it with
    # stderr closed keeps bash from printing its own crash message.
    ( cd "$wd" && exec perl -e 'alarm shift @ARGV; exec @ARGV or die "cannot run $ARGV[0]: $!\n"' \
        "$TOOL_TIMEOUT" "$BIN/$name" < /dev/null > _stdout.txt 2> _stderr.txt ) &
    { wait $!; } 2> /dev/null || rc=$?
    t1=$(now)
    [ -s "$wd/_stderr.txt" ] || rm -f "$wd/_stderr.txt"
    note=ok
    if [ "$rc" -eq 142 ]; then note="killed: timeout ${TOOL_TIMEOUT}s"
    elif [ "$rc" -ge 128 ]; then note="signal $((rc - 128))"
    fi
    printf '%s\t%s\t%s\t%s\n' "$name" "$rc" "$(perl -e "printf '%.2f', $t1 - $t0")" "$note" >> "$status"
    echo "  $name exit $rc"
done
