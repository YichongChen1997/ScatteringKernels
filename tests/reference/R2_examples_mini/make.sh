#!/usr/bin/env bash
# R2 reference answers: examples/Kerogen (methane between two kerogen
# surfaces, 423 K, Poiseuille flow; the case of the Fuel 2022 paper), mini size.
#
# Regenerates everything from a checkout of tag v2.0.0:
#   init   : examples/Kerogen/initialisation/combine.cpp (unchanged) -> data.dat
#   md     : edited copies of examples/Kerogen/in.equil and in.meas, mpirun -np $NP lmp
#   bin    : all tools/*.cpp compiled with plain "g++ -o name name.cpp" (no flags)
#   tools  : every tool on the first $NMINI frames of the gas dump
#   repeat : every tool again, in a fresh directory; stop unless byte-identical
#   collect: copy the small files into OUT (default: this script's directory)
#
# Usage: make.sh V2_CHECKOUT [RUN_DIR] [OUT_DIR]
#   V2_CHECKOUT  clean checkout of tag v2.0.0 (read only)
#   RUN_DIR      where big outputs go, about 1.6 GB (default: $TMPDIR/sk_R2_run)
#   OUT_DIR      where the stored reference files go (default: script directory)
# Environment: NP (default 4), LMP (default lmp), MPIRUN (default mpirun),
#   CXX (default g++), STAGES (default "init md bin tools repeat collect"),
#   TOOL_TIMEOUT seconds per tool (default 600).
set -euo pipefail
export LC_ALL=C
export OMP_NUM_THREADS=1

usage() { sed -n '2,19p' "$0" >&2; exit 2; }
[ $# -ge 1 ] || usage

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
V2=$(cd "$1" && pwd)
RUN=${2:-${TMPDIR:-/tmp}/sk_R2_run}
OUT=${3:-$HERE}
mkdir -p "$RUN" "$OUT"
RUN=$(cd "$RUN" && pwd)
OUT=$(cd "$OUT" && pwd)

NP=${NP:-4}
LMP=${LMP:-lmp}
MPIRUN=${MPIRUN:-mpirun}
CXX=${CXX:-g++}
STAGES=${STAGES:-"init md bin tools repeat collect"}
TOOL_TIMEOUT=${TOOL_TIMEOUT:-600}

NMINI=201                        # frames kept in the regression dump (all of them here)
MINI_GZ_LIMIT=1500000            # bytes; gzip'd regression dump must be smaller
STORE_LIMIT=2000000              # bytes; larger tool outputs are stored as sha256 only
RSS_LIMIT_KB=$((8 * 1024 * 1024))  # kill a tool whose resident memory passes 8 GB

EXDIR=$V2/examples/Kerogen
KEROGEN_XYZ=EFK_50A_0.80.xyz     # the structure combine.cpp opens by name
BIN=$RUN/bin
GZNAME=dump_meas_gas.first$NMINI.lammpstrj.gz
GZ=$RUN/$GZNAME                 # the regression dump, compressed (made by stage tools)

log() { printf '[make.sh %s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
die() { printf '[make.sh] ERROR: %s\n' "$*" >&2; exit 1; }
now() { perl -MTime::HiRes=time -e 'printf "%.3f\n", time'; }
sha() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
sha_stdin() { sha - | awk '{print $1}'; }
fsize() { wc -c < "$1" | tr -d ' '; }
has_stage() { case " $STAGES " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }
nframes() { grep -c '^ITEM: TIMESTEP' "$1"; }

# edit_file SRC DST  LINE OLD NEW  [LINE OLD NEW ...]
# Replace whole lines by number, but only if the line is exactly OLD.
# Everything else (including a missing final newline) is kept byte for byte.
edit_file() {
    local src=$1 dst=$2 tmp
    shift 2
    tmp="$dst.tmp.$$"
    cp "$src" "$tmp"
    while [ $# -gt 0 ]; do
        EDIT_LINE=$1 EDIT_OLD=$2 EDIT_NEW=$3 perl -pe '
            if ($. == $ENV{EDIT_LINE}) {
                my $nl = s/(\r?\n)\z// ? $1 : "";
                die "edit_file: line $. is <$_>, expected <$ENV{EDIT_OLD}>\n"
                    unless $_ eq $ENV{EDIT_OLD};
                $_ = $ENV{EDIT_NEW} . $nl;
            }' "$tmp" > "$tmp.2"
        mv "$tmp.2" "$tmp"
        shift 3
    done
    mv "$tmp" "$dst"
}

# Specification.dat read by every tools/pp_meas_*.cpp (5 lines). Values for
# the Kerogen decks: time step 0.5 fs, gas dump every 500 steps, virtual plane
# at the 15 A gas-wall cut-off, channel height H = 100 A (combine.cpp),
# gas and wall at 423 K, methane 16.043 g/mol, kerogen carbon 12.0107 g/mol.
write_spec() {
    printf '%d                  # nTimeSteps (frames to read)\n' "$1"
    printf '0.5  500  0          # deltaT[fs]  tSkip(dump every)  skipTimeStep(frames)\n'
    printf '15   100             # rCut (virtual plane z, A)  H [A]\n'
    printf '423  423             # Tg  Tw [K]\n'
    printf '16.043  12.0107      # mG  mW [g/mol]\n'
}

# /usr/bin/time flavour: BSD/macOS "-l" or GNU "-v"; none -> no rusage report.
TIME_FLAG=""
if /usr/bin/time -l true > /dev/null 2>&1; then TIME_FLAG=-l
elif /usr/bin/time -v true > /dev/null 2>&1; then TIME_FLAG=-v
fi

# run_tool NAME WORKDIR STATUSFILE -> appends one line to STATUSFILE
# Each tool runs alone in WORKDIR (which holds links to the gas dump and
# Specification.dat); the tool's stdout -> _stdout.txt, its stderr ->
# _stderr.txt (removed if empty); /usr/bin/time report -> _time.txt and
# time's own messages -> _time_stderr.txt (both are not reference data).
run_tool() {
    local name=$1 wd=$2 status=$3
    local t0 t1 pid child rss peak=0 rc=0 note="" elapsed
    t0=$(now)
    if [ -n "$TIME_FLAG" ]; then
        ( cd "$wd" && exec /usr/bin/time "$TIME_FLAG" -o _time.txt \
            sh -c 'exec "$0" < /dev/null > _stdout.txt 2> _stderr.txt' "$BIN/$name" 2> _time_stderr.txt ) &
    else
        ( cd "$wd" && exec "$BIN/$name" < /dev/null > _stdout.txt 2> _stderr.txt ) &
    fi
    pid=$!
    while kill -0 "$pid" 2>/dev/null; do
        child=$(pgrep -P "$pid" 2>/dev/null | head -1 || true)
        [ -n "$child" ] || child=$pid
        if [ -n "$child" ]; then
            rss=$(ps -o rss= -p "$child" 2>/dev/null | tr -d ' ' || true)
            if [ -n "$rss" ] && [ "$rss" -gt "$peak" ]; then peak=$rss; fi
            if [ -n "$rss" ] && [ "$rss" -gt "$RSS_LIMIT_KB" ]; then
                note="killed: resident memory > 8 GB"; kill -9 "$child" 2>/dev/null || true
            fi
            elapsed=$(perl -e "printf '%d', $(now) - $t0")
            if [ "$elapsed" -gt "$TOOL_TIMEOUT" ]; then
                note="killed: timeout ${TOOL_TIMEOUT}s"; kill -9 "$child" 2>/dev/null || true
            fi
        fi
        sleep 0.2
    done
    wait "$pid" || rc=$?
    t1=$(now)
    [ -s "$wd/_stderr.txt" ] || rm -f "$wd/_stderr.txt"
    # exact peak from getrusage when available (bytes on macOS, kB with GNU time)
    local maxrss=""
    if [ -f "$wd/_time.txt" ]; then
        maxrss=$(awk '/maximum resident set size/ {print $1; exit}
                      /Maximum resident set size/ {print $NF * 1024; exit}' "$wd/_time.txt")
    fi
    [ -n "$maxrss" ] || maxrss=$((peak * 1024))
    if [ -z "$note" ] && [ "$rc" -ge 128 ]; then note="signal $((rc - 128))"; fi
    printf '%s\t%s\t%.2f\t%.0f\t%s\n' "$name" "$rc" "$(perl -e "print $t1 - $t0")" \
        "$(perl -e "print $maxrss / 1048576")" "${note:-ok}" >> "$status"
    log "  $name rc=$rc ${note:-}"
}

# run_all_tools DIR : DIR holds dump_meas_gas.lammpstrj and Specification.dat.
# Each tool gets its own subdirectory DIR/tools/<name>. The Kerogen decks
# write no wall dump, so pp_meas_wallTemp finds no dump_meas_wall.lammpstrj.
run_all_tools() {
    local dir=$1 src name wd
    local status=$dir/tool_status.tsv
    printf 'tool\texit_code\twall_s\tmax_rss_MB\tnote\n' > "$status"
    rm -rf "$dir/tools"
    for src in "$V2"/tools/*.cpp; do
        name=$(basename "$src" .cpp)
        wd=$dir/tools/$name
        mkdir -p "$wd"
        ln -s ../../dump_meas_gas.lammpstrj "$wd/dump_meas_gas.lammpstrj"
        ln -s ../../Specification.dat "$wd/Specification.dat"
        if [ -x "$BIN/$name" ]; then
            run_tool "$name" "$wd" "$status"
        else
            printf '%s\tNA\tNA\tNA\tnot built\n' "$name" >> "$status"
        fi
    done
    # stdout tails, for the record
    : > "$dir/tool_stdout_tail.txt"
    for wd in "$dir"/tools/*; do
        {
            printf '==== %s (last 4 stdout lines) ====\n' "$(basename "$wd")"
            tail -n 4 "$wd/_stdout.txt" 2>/dev/null || true
            if [ -s "$wd/_stderr.txt" ]; then
                printf -- '---- stderr (last 4 lines) ----\n'; tail -n 4 "$wd/_stderr.txt"
            fi
        } >> "$dir/tool_stdout_tail.txt"
    done
}

# sha256 of every tool output (and its stdout) under DIR/tools, relative paths,
# skipping the input links and the timing report.
output_sums() {
    local dir=$1
    ( cd "$dir" && find tools -type f ! -name _time.txt ! -name _time_stderr.txt | sort | while read -r f; do sha "$f"; done )
}

# prepare_ana DIR : DIR/run gets the regression dump (from the .gz, exactly
# as a user of the stored set would have it) and Specification.dat.
prepare_ana() {
    local dir=$1
    rm -rf "$dir"; mkdir -p "$dir/run"
    gunzip -c "$GZ" > "$dir/run/dump_meas_gas.lammpstrj"
    write_spec "$NMINI" > "$dir/run/Specification.dat"
}

########################################################################
log "V2=$V2 RUN=$RUN OUT=$OUT NP=$NP STAGES=$STAGES"

if has_stage init; then
    log "stage init"
    rm -rf "$RUN/init"; mkdir -p "$RUN/init"; cd "$RUN/init"
    cp "$EXDIR/initialisation/combine.cpp" "$EXDIR/initialisation/$KEROGEN_XYZ" .
    "$CXX" -o combine combine.cpp > combine.build.log 2>&1 || { cat combine.build.log >&2; die "combine.cpp did not compile"; }
    ./combine > combine.stdout
    [ -s data.dat ] || die "combine wrote no data.dat"
    grep -Eq '^6189[[:space:]]+atoms$' data.dat || die "data.dat does not have 6189 atoms (see combine.stdout)"
fi

if has_stage md; then
    log "stage md (mpirun -np $NP)"
    rm -rf "$RUN/md"; mkdir -p "$RUN/md"; cd "$RUN/md"
    cp "$RUN/init/data.dat" .
    edit_file "$EXDIR/in.equil" in.equil \
        17 'processors      16 16 2' 'processors      * * *' \
        74 'run             200000 every 100000 "write_restart restart.*"' \
           'run             20000 every 10000 "write_restart restart.*"'
    edit_file "$EXDIR/in.meas" in.meas \
        17 'processors      16 16 2 ' 'processors      * * * ' \
        90 'run             10000000 every 1000000 "write_restart restart.*"' \
           'run             100000 every 10000 "write_restart restart.*"'
    t0=$(now)
    "$MPIRUN" -np "$NP" "$LMP" -in in.equil > screen.equil 2>&1 || { tail -20 screen.equil >&2; die "in.equil failed"; }
    t1=$(now)
    "$MPIRUN" -np "$NP" "$LMP" -in in.meas > screen.meas 2>&1 || { tail -20 screen.meas >&2; die "in.meas failed"; }
    t2=$(now)
    printf 'equil_wall_s %.1f\nmeas_wall_s %.1f\n' "$(perl -e "print $t1-$t0")" "$(perl -e "print $t2-$t1")" > md_times.txt
    cat md_times.txt
    for f in log.equil log.meas; do
        grep -q '^Total wall time' "$f" || die "$f has no 'Total wall time' line"
    done
    n=$(nframes dump_meas_gas.lammpstrj)
    [ "$n" -ge "$NMINI" ] || die "gas dump has $n frames, need at least $NMINI"
    log "gas dump: $n frames, $(fsize dump_meas_gas.lammpstrj) bytes"
fi

if has_stage bin; then
    log "stage bin ($CXX -o name name.cpp, no flags)"
    rm -rf "$BIN"; mkdir -p "$BIN"
    : > "$BIN/build_status.tsv"
    for src in "$V2"/tools/*.cpp; do
        name=$(basename "$src" .cpp)
        rc=0
        ( cd "$BIN" && "$CXX" -o "$name" "$src" ) > "$BIN/$name.build.log" 2>&1 || rc=$?
        printf '%s\t%s\t%s warning lines\n' "$name" "$rc" "$(grep -c 'warning:' "$BIN/$name.build.log" || true)" >> "$BIN/build_status.tsv"
    done
    cat "$BIN/build_status.tsv"
fi

if has_stage tools; then
    log "stage tools (first $NMINI frames)"
    M=$RUN/ana
    awk -v N="$NMINI" '/^ITEM: TIMESTEP/{n++} n>N{exit} {print}' \
        "$RUN/md/dump_meas_gas.lammpstrj" | gzip -9n > "$GZ"
    gz=$(fsize "$GZ")
    [ "$gz" -lt "$MINI_GZ_LIMIT" ] || die "regression dump gz is $gz bytes >= $MINI_GZ_LIMIT"
    prepare_ana "$M"
    n=$(nframes "$M/run/dump_meas_gas.lammpstrj")
    [ "$n" -eq "$NMINI" ] || die "regression dump has $n frames, expected $NMINI"
    run_all_tools "$M/run"
    output_sums "$M/run" > "$M/run/OUTPUT_SHA256SUMS"
fi

if has_stage repeat; then
    log "stage repeat (all tools again in a fresh directory)"
    R=$RUN/ana_repeat
    prepare_ana "$R"
    run_all_tools "$R/run"
    output_sums "$R/run" > "$R/run/OUTPUT_SHA256SUMS"
    if cmp -s "$RUN/ana/run/OUTPUT_SHA256SUMS" "$R/run/OUTPUT_SHA256SUMS"; then
        printf 'identical: %s files\n' "$(wc -l < "$R/run/OUTPUT_SHA256SUMS" | tr -d ' ')" > "$R/REPEAT_RESULT"
        log "repeat: $(cat "$R/REPEAT_RESULT")"
    else
        diff "$RUN/ana/run/OUTPUT_SHA256SUMS" "$R/run/OUTPUT_SHA256SUMS" >&2 || true
        die "tool outputs differ between two runs on the same dump"
    fi
fi

if has_stage collect; then
    log "stage collect -> $OUT"
    M=$RUN/ana
    # inputs (the edited copies actually run) and their diffs against v2.0.0
    rm -rf "$OUT/inputs"; mkdir -p "$OUT/inputs"
    cp "$RUN/md/in.equil" "$RUN/md/in.meas" "$OUT/inputs/"
    {
        for f in in.equil in.meas; do
            printf '### diff v2.0.0:examples/Kerogen/%s inputs/%s\n' "$f" "$f"
            diff "$EXDIR/$f" "$OUT/inputs/$f" || true
        done
    } > "$OUT/inputs/EDITS.diff"

    # the full run: hashes only (paths relative to RUN)
    ( cd "$RUN" && sha init/data.dat md/dump_equil_gas.lammpstrj md/dump_meas_gas.lammpstrj ) > "$OUT/RUN_SHA256SUMS"

    # mini regression set: inputs, every output (<2 MB) or its hash
    rm -rf "$OUT/mini"; mkdir -p "$OUT/mini/outputs"
    cp "$GZ" "$M/run/Specification.dat" "$M/run/tool_status.tsv" "$M/run/tool_stdout_tail.txt" "$OUT/mini/"
    # SHA256SUMS lists only files stored here, so that
    #   (cd mini && shasum -a 256 -c SHA256SUMS)
    # passes as it is. Outputs too large to store keep their sha256 and size
    # in outputs/NOT_STORED.txt. The .gz bytes depend on the gzip build;
    # DUMP_CONTENT.sha256 has the sha256 of the uncompressed dump.
    printf '# sha256  bytes  path (relative to outputs/), not stored because > %s bytes\n' "$STORE_LIMIT" \
        > "$OUT/mini/outputs/NOT_STORED.txt"
    : > "$OUT/mini/stored_outputs.sha256"
    ( cd "$M/run/tools" && find . -type f ! -name _time.txt ! -name _time_stderr.txt | sort ) | while read -r f; do
        f=${f#./}
        line=$(grep -F "  tools/$f" "$M/run/OUTPUT_SHA256SUMS" | awk -v p="tools/$f" '$2 == p')
        [ -n "$line" ] || die "no sha256 for tools/$f"
        if [ "$(fsize "$M/run/tools/$f")" -lt "$STORE_LIMIT" ]; then
            mkdir -p "$OUT/mini/outputs/$(dirname "$f")"
            cp "$M/run/tools/$f" "$OUT/mini/outputs/$f"
            printf '%s  outputs/%s\n' "${line%%  *}" "$f" >> "$OUT/mini/stored_outputs.sha256"
        else
            printf '%s  %s  %s\n' "${line%%  *}" "$(fsize "$M/run/tools/$f")" "$f" >> "$OUT/mini/outputs/NOT_STORED.txt"
        fi
    done
    {
        ( cd "$OUT/mini" && sha "$GZNAME" Specification.dat )
        cat "$OUT/mini/stored_outputs.sha256"
    } > "$OUT/mini/SHA256SUMS"
    rm -f "$OUT/mini/stored_outputs.sha256"
    printf '%s  dump_meas_gas.first%s.lammpstrj\n' \
        "$(gzip -dc "$OUT/mini/$GZNAME" | sha_stdin)" "$NMINI" \
        > "$OUT/mini/DUMP_CONTENT.sha256"

    # facts about the run; the tool outputs themselves are not physical numbers (KI-8)
    D=$RUN/md/dump_meas_gas.lammpstrj
    natoms=$(awk 'prev ~ /^ITEM: NUMBER OF ATOMS/ {print; exit} {prev = $0}' "$D")
    cols=$(awk '/^ITEM: ATOMS/ {print NF - 2; exit}' "$D")
    first=$(awk 'prev ~ /^ITEM: TIMESTEP/ {print; exit} {prev = $0}' "$D")
    last=$(awk 'prev ~ /^ITEM: TIMESTEP/ {s = $0} {prev = $0} END {print s}' "$D")
    total=$(awk '$2 == "atoms" && NF == 2 {print $1; exit}' "$RUN/init/data.dat")
    for v in "$natoms" "$cols" "$first" "$last" "$total"; do
        case "$v" in ''|*[!0-9]*) die "could not read a number from the run files: <$v>" ;; esac
    done
    {
        printf '{\n'
        printf '  "case": "R2_examples_mini (examples/Kerogen)",\n'
        printf '  "source": "tag v2.0.0 (%s)",\n' "$(git -C "$V2" rev-parse HEAD 2>/dev/null || echo unknown)"
        printf '  "np": %s,\n' "$NP"
        printf '  "atoms_total": %s,\n' "$total"
        printf '  "gas_atoms_per_frame": %s,\n' "$natoms"
        printf '  "dump_columns": %s,\n' "$cols"
        printf '  "columns_the_tools_read": 16,\n'
        printf '  "frames_full_dump": %s,\n' "$(nframes "$D")"
        printf '  "frames_regression_dump": %s,\n' "$NMINI"
        printf '  "first_and_last_timestep": [%s, %s],\n' "$first" "$last"
        printf '  "bytes": {"data.dat": %s, "dump_equil_gas.lammpstrj": %s, "dump_meas_gas.lammpstrj": %s},\n' \
            "$(fsize "$RUN/init/data.dat")" "$(fsize "$RUN/md/dump_equil_gas.lammpstrj")" "$(fsize "$D")"
        printf '  "note": "the dump has 11 columns and every tool reads 16 (KI-8), so tool outputs are regression data only, not physical results"\n'
        printf '}\n'
    } > "$OUT/mini/expected.json"

    # provenance (changes on every run: date, host load)
    {
        printf 'date: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        printf 'uname: %s\n' "$(uname -srm)"
        printf 'v2 checkout: %s\n' "$(git -C "$V2" describe --tags --always 2>/dev/null || echo unknown) $(git -C "$V2" rev-parse HEAD 2>/dev/null || true)"
        printf 'np: %s\nOMP_NUM_THREADS: %s\n' "$NP" "$OMP_NUM_THREADS"
        printf -- '--- %s -h | head -3\n' "$(basename "$LMP")"; "$LMP" -h | head -3
        printf -- '--- %s --version\n' "$MPIRUN"; "$MPIRUN" --version 2>&1 | head -1
        printf -- '--- %s --version\n' "$CXX"; "$CXX" --version 2>&1 | head -2
        printf -- '--- processor grid\n'; grep -h 'MPI processor grid' "$RUN/md/log.equil" "$RUN/md/log.meas" 2>/dev/null || true
        printf -- '--- md wall times\n'; cat "$RUN/md/md_times.txt" 2>/dev/null || true
        printf -- '--- tool build status (name rc warnings)\n'; cat "$BIN/build_status.tsv" 2>/dev/null || true
        printf -- '--- repeat stage\n'; cat "$RUN/ana_repeat/REPEAT_RESULT" 2>/dev/null || echo 'not run'
        printf -- '--- v2.0.0 sources used (sha256)\n'
        ( cd "$V2" && sha examples/Kerogen/initialisation/combine.cpp \
            "examples/Kerogen/initialisation/$KEROGEN_XYZ" \
            examples/Kerogen/in.equil examples/Kerogen/in.meas tools/*.cpp )
    } > "$OUT/provenance.txt"
fi

log "done"
