#!/usr/bin/env bash
# R3 reference answers: thesis Chapter-3 channel case (explicit layered Pt walls,
# Ar between parallel walls, 300 K, equilibrium MD), mini size.
#
# Regenerates everything from a checkout of tag v2.0.0:
#   init  : combine.cpp + edited Parameters.dat -> data.dat
#   md    : edited in.equil / in.meas, mpirun -np $NP lmp
#   bin   : all tools/*.cpp compiled with plain "g++ -o name name.cpp" (no flags)
#   full  : every tool on the full gas dump (all frames)
#   mini  : every tool on the first $NMINI frames (the byte-exact regression set)
#   collect: copy the small files into OUT (default: this script's directory)
#
# Usage: make.sh V2_CHECKOUT [RUN_DIR] [OUT_DIR]
#   V2_CHECKOUT  clean checkout of tag v2.0.0 (read only)
#   RUN_DIR      where big outputs go, about 2 GB (default: $TMPDIR/sk_R3_run)
#   OUT_DIR      where the stored reference files go (default: script directory)
# Environment: NP (default 4), LMP (default lmp), MPIRUN (default mpirun),
#   CXX (default g++), STAGES (default "init md bin full mini collect"),
#   TOOL_TIMEOUT seconds per tool (default 3600).
set -euo pipefail
export LC_ALL=C
export OMP_NUM_THREADS=1

usage() { sed -n '2,19p' "$0" >&2; exit 2; }
[ $# -ge 1 ] || usage

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
V2=$(cd "$1" && pwd)
RUN=${2:-${TMPDIR:-/tmp}/sk_R3_run}
OUT=${3:-$HERE}
mkdir -p "$RUN" "$OUT"
RUN=$(cd "$RUN" && pwd)
OUT=$(cd "$OUT" && pwd)

NP=${NP:-4}
LMP=${LMP:-lmp}
MPIRUN=${MPIRUN:-mpirun}
CXX=${CXX:-g++}
STAGES=${STAGES:-"init md bin full mini collect"}
TOOL_TIMEOUT=${TOOL_TIMEOUT:-3600}

NMINI=80                         # frames kept in the regression dump
MINI_GZ_LIMIT=1500000            # bytes; gzip'd truncated dump must be smaller
STORE_LIMIT=2000000              # bytes; larger tool outputs are stored as sha256 only
RSS_LIMIT_KB=$((8 * 1024 * 1024))  # kill a tool whose resident memory passes 8 GB

WALLDIR=$V2/initialisation/TypeOfWalls/Explicit_layered
DECKDIR=$V2/LAMMPS_input/Explicit_layered
BIN=$RUN/bin

log() { printf '[make.sh %s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
now() { perl -MTime::HiRes=time -e 'printf "%.3f\n", time'; }
sha() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
sha_stdin() { sha - | awk '{print $1}'; }
fsize() { wc -c < "$1" | tr -d ' '; }
has_stage() { case " $STAGES " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

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

# Specification.dat read by every tools/pp_meas_*.cpp (5 lines).
write_spec() {
    printf '%d                  # nTimeSteps (frames to read)\n' "$1"
    printf '1.0  200  0          # deltaT[fs]  tSkip(dump every)  skipTimeStep(frames)\n'
    printf '12   30              # rCut (virtual plane z, A)  H [A]\n'
    printf '300  300             # Tg  Tw [K]\n'
    printf '39.948  195.084      # mG  mW [g/mol]\n'
}

# /usr/bin/time flavour: BSD/macOS "-l" or GNU "-v"; none -> no rusage report.
TIME_FLAG=""
if /usr/bin/time -l true > /dev/null 2>&1; then TIME_FLAG=-l
elif /usr/bin/time -v true > /dev/null 2>&1; then TIME_FLAG=-v
fi

# run_tool NAME WORKDIR STATUSFILE -> appends one line to STATUSFILE
# Each tool runs alone in WORKDIR (which holds links to the dumps and
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

# run_all_tools DIR : DIR holds dump_meas_gas.lammpstrj, dump_meas_wall.lammpstrj,
# Specification.dat. Each tool gets its own subdirectory DIR/tools/<name>.
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
        ln -s ../../dump_meas_wall.lammpstrj "$wd/dump_meas_wall.lammpstrj"
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

# ACs numbers: row 1 of Alpha*.txt; column 4 = general formula (thesis eq. 3.6,
# all collisions), column 5 = least squares (eq. 3.7, all collisions).
# Columns 4 and 5 are the same on every row; check that.
ac_json() {
    local dir=$1 f key sep=""
    printf '{'
    for key in TMAC:AlphaTx NMAC:AlphaN NEAC:AlphaEn EAC:AlphaE TEAC_x:AlphaEx; do
        f=$dir/tools/pp_meas_ACs/${key#*:}.txt
        if [ "$(awk '{print $4, $5}' "$f" | sort -u | wc -l | tr -d ' ')" != 1 ]; then
            echo "ac_json: columns 4-5 of $f differ between rows" >&2; exit 1
        fi
        printf '%s"%s": {"general": %s, "least_squares": %s}' "$sep" "${key%%:*}" \
            "$(awk 'NR==1{print $4}' "$f")" "$(awk 'NR==1{print $5}' "$f")"
        sep=", "
    done
    printf '}'
}

collisions_json() {
    local o=$1/tools/pp_meas_ACs/_stdout.txt
    printf '{"started": %s, "ended": %s, "bottom": %s, "top": %s}' \
        "$(awk -F': ' '/Number of Collisions Started/{print $2}' "$o")" \
        "$(awk -F': ' '/Number of Collisions Ended/{print $2}' "$o")" \
        "$(awk -F': ' '/No. of Collisions at bottom/{print $2}' "$o")" \
        "$(awk -F': ' '/No. of Collisions at top/{print $2}' "$o")"
}

# A value that awk failed to find leaves an empty field, e.g. '"ended": ,'.
# Stop instead of writing broken JSON.
require_json_values() {
    case "$1" in
        *': ,'* | *': }'* | *':,'* | *':}'* | '')
            echo "missing value in: $1" >&2; exit 1 ;;
    esac
}

########################################################################
log "V2=$V2 RUN=$RUN OUT=$OUT NP=$NP STAGES=$STAGES"

if has_stage init; then
    log "stage init"
    rm -rf "$RUN/init"; mkdir -p "$RUN/init"; cd "$RUN/init"
    cp "$WALLDIR/combine.cpp" "$WALLDIR/unitCell.dat" .
    edit_file "$WALLDIR/Parameters.dat" Parameters.dat \
        1 '0    200          # xLo,  xHi'  '0    40           # xLo,  xHi' \
        2 '0    200          # yLo,  yHi ' '0    40           # yLo,  yHi '
    "$CXX" -o combine combine.cpp > combine.build.log 2>&1
    ./combine > combine.stdout
fi

if has_stage md; then
    log "stage md (mpirun -np $NP)"
    rm -rf "$RUN/md"; mkdir -p "$RUN/md"; cd "$RUN/md"
    cp "$RUN/init/data.dat" .
    edit_file "$DECKDIR/in.equil" in.equil \
        16  'processors      8 8 2' 'processors      * * *' \
        106 'run             100000 every 100000 "write_restart restart.*"' \
            'run             50000 every 50000 "write_restart restart.*"'
    edit_file "$DECKDIR/in.meas" in.meas \
        16  'processors      8 8 2' 'processors      * * *' \
        51  'pair_coeff   1 1 ${epsilonGas}   ${sigmaGas}   ${rCut}' \
            '#pair_coeff   1 1 ${epsilonGas}   ${sigmaGas}   ${rCut}' \
        52  '#pair_coeff   1 1 0.0000 ${sigmaGas} 0.1' \
            'pair_coeff   1 1 0.0000 ${sigmaGas} 0.1' \
        106 'run             6000000 every 6000000 "write_restart restart.*"' \
            'run             230000 every 230000 "write_restart restart.*"'
    t0=$(now)
    "$MPIRUN" -np "$NP" "$LMP" -in in.equil > screen.equil 2>&1
    t1=$(now)
    "$MPIRUN" -np "$NP" "$LMP" -in in.meas > screen.meas 2>&1
    t2=$(now)
    printf 'equil_wall_s %.1f\nmeas_wall_s %.1f\n' "$(perl -e "print $t1-$t0")" "$(perl -e "print $t2-$t1")" > md_times.txt
    cat md_times.txt
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

if has_stage full; then
    log "stage full (all frames)"
    A=$RUN/ana_full
    rm -rf "$A"; mkdir -p "$A"
    ln -s ../md/dump_meas_gas.lammpstrj "$A/dump_meas_gas.lammpstrj"
    ln -s ../md/dump_meas_wall.lammpstrj "$A/dump_meas_wall.lammpstrj"
    NFULL=$(grep -c '^ITEM: TIMESTEP' "$RUN/md/dump_meas_gas.lammpstrj")
    write_spec "$NFULL" > "$A/Specification.dat"
    run_all_tools "$A"
    output_sums "$A" > "$A/OUTPUT_SHA256SUMS"
fi

if has_stage mini; then
    log "stage mini (first $NMINI frames)"
    M=$RUN/ana_mini
    rm -rf "$M"; mkdir -p "$M/run"
    awk -v N="$NMINI" '/^ITEM: TIMESTEP/{n++} n>N{exit} {print}' \
        "$RUN/md/dump_meas_gas.lammpstrj" | gzip -9n > "$M/dump_meas_gas.first$NMINI.lammpstrj.gz"
    gzip -9n < "$RUN/md/dump_meas_wall.lammpstrj" > "$M/dump_meas_wall.lammpstrj.gz"
    gz=$(fsize "$M/dump_meas_gas.first$NMINI.lammpstrj.gz")
    [ "$gz" -lt "$MINI_GZ_LIMIT" ] || { echo "truncated dump gz is $gz bytes >= $MINI_GZ_LIMIT" >&2; exit 1; }
    # run from the compressed files, exactly as a user of the stored set would
    gunzip -c "$M/dump_meas_gas.first$NMINI.lammpstrj.gz" > "$M/run/dump_meas_gas.lammpstrj"
    gunzip -c "$M/dump_meas_wall.lammpstrj.gz" > "$M/run/dump_meas_wall.lammpstrj"
    write_spec "$NMINI" > "$M/run/Specification.dat"
    run_all_tools "$M/run"
    output_sums "$M/run" > "$M/run/OUTPUT_SHA256SUMS"
fi

if has_stage collect; then
    log "stage collect -> $OUT"
    A=$RUN/ana_full; M=$RUN/ana_mini
    # inputs (the edited copies actually run) and their diffs against v2.0.0
    rm -rf "$OUT/inputs"; mkdir -p "$OUT/inputs"
    cp "$RUN/init/Parameters.dat" "$RUN/md/in.equil" "$RUN/md/in.meas" "$OUT/inputs/"
    {
        for p in "initialisation/TypeOfWalls/Explicit_layered/Parameters.dat:Parameters.dat" \
                 "LAMMPS_input/Explicit_layered/in.equil:in.equil" \
                 "LAMMPS_input/Explicit_layered/in.meas:in.meas"; do
            printf '### diff v2.0.0:%s inputs/%s\n' "${p%%:*}" "${p#*:}"
            diff "$V2/${p%%:*}" "$OUT/inputs/${p#*:}" || true
        done
    } > "$OUT/inputs/EDITS.diff"

    # full run: hashes and numbers only
    rm -rf "$OUT/full"; mkdir -p "$OUT/full"
    cp "$A/Specification.dat" "$A/tool_status.tsv" "$A/tool_stdout_tail.txt" "$A/OUTPUT_SHA256SUMS" "$OUT/full/"
    ( cd "$RUN" && sha init/data.dat md/dump_meas_gas.lammpstrj md/dump_meas_wall.lammpstrj ) > "$OUT/full/SHA256SUMS"
    NFULL=$(awk 'NR==1{print $1}' "$A/Specification.dat")
    # assigned first so that a failure stops the script (set -e)
    col_full=$(collisions_json "$A"); require_json_values "$col_full"
    ac_full=$(ac_json "$A"); require_json_values "$ac_full"
    {
        printf '{\n'
        printf '  "case": "R3_channel_mini",\n'
        printf '  "source": "tag v2.0.0 (%s)",\n' "$(git -C "$V2" rev-parse HEAD 2>/dev/null || echo unknown)"
        printf '  "np": %s,\n' "$NP"
        printf '  "frames": %s,\n' "$NFULL"
        printf '  "sha256": {"data.dat": "%s", "dump_meas_gas.lammpstrj": "%s", "dump_meas_wall.lammpstrj": "%s"},\n' \
            $(awk '{print $1}' "$OUT/full/SHA256SUMS")
        printf '  "dump_bytes": {"dump_meas_gas.lammpstrj": %s, "dump_meas_wall.lammpstrj": %s},\n' \
            "$(fsize "$RUN/md/dump_meas_gas.lammpstrj")" "$(fsize "$RUN/md/dump_meas_wall.lammpstrj")"
        printf '  "collisions": %s,\n' "$col_full"
        printf '  "accommodation_coefficients_bottom_wall": %s,\n' "$ac_full"
        printf '  "ac_columns": "Alpha*.txt row 1: column 4 = general formula (eq. 3.6), column 5 = least squares (eq. 3.7), both over all collisions",\n'
        printf '  "thesis_fig_3_4": {"TMAC": 0.49, "NMAC": 0.67, "NEAC": 0.63, "EAC": 0.26}\n'
        printf '}\n'
    } > "$OUT/full/expected.json"

    # mini regression set: inputs, every output (<2 MB) or its hash
    rm -rf "$OUT/mini"; mkdir -p "$OUT/mini/outputs"
    cp "$M/dump_meas_gas.first$NMINI.lammpstrj.gz" "$M/dump_meas_wall.lammpstrj.gz" \
       "$M/run/Specification.dat" "$M/run/tool_status.tsv" "$M/run/tool_stdout_tail.txt" "$OUT/mini/"
    # SHA256SUMS lists only files stored here, so that
    #   (cd mini && shasum -a 256 -c SHA256SUMS)
    # passes as it is. Outputs too large to store keep their sha256 and size
    # in outputs/NOT_STORED.txt. The .gz bytes depend on the gzip build;
    # DUMP_CONTENT.sha256 has the sha256 of the uncompressed dumps.
    printf '# sha256  bytes  path (relative to outputs/), not stored because > %s bytes\n' "$STORE_LIMIT" \
        > "$OUT/mini/outputs/NOT_STORED.txt"
    : > "$OUT/mini/stored_outputs.sha256"
    ( cd "$M/run/tools" && find . -type f ! -name _time.txt ! -name _time_stderr.txt | sort ) | while read -r f; do
        f=${f#./}
        line=$(grep -F "  tools/$f" "$M/run/OUTPUT_SHA256SUMS" | awk -v p="tools/$f" '$2 == p')
        [ -n "$line" ] || { echo "no sha256 for tools/$f" >&2; exit 1; }
        if [ "$(fsize "$M/run/tools/$f")" -lt "$STORE_LIMIT" ]; then
            mkdir -p "$OUT/mini/outputs/$(dirname "$f")"
            cp "$M/run/tools/$f" "$OUT/mini/outputs/$f"
            printf '%s  outputs/%s\n' "${line%%  *}" "$f" >> "$OUT/mini/stored_outputs.sha256"
        else
            printf '%s  %s  %s\n' "${line%%  *}" "$(fsize "$M/run/tools/$f")" "$f" >> "$OUT/mini/outputs/NOT_STORED.txt"
        fi
    done
    {
        ( cd "$OUT/mini" && sha "dump_meas_gas.first$NMINI.lammpstrj.gz" dump_meas_wall.lammpstrj.gz Specification.dat )
        cat "$OUT/mini/stored_outputs.sha256"
    } > "$OUT/mini/SHA256SUMS"
    rm -f "$OUT/mini/stored_outputs.sha256"
    {
        printf '%s  dump_meas_gas.first%s.lammpstrj\n' "$(gzip -dc "$OUT/mini/dump_meas_gas.first$NMINI.lammpstrj.gz" | sha_stdin)" "$NMINI"
        printf '%s  dump_meas_wall.lammpstrj\n' "$(gzip -dc "$OUT/mini/dump_meas_wall.lammpstrj.gz" | sha_stdin)"
    } > "$OUT/mini/DUMP_CONTENT.sha256"
    col_mini=$(collisions_json "$M/run"); require_json_values "$col_mini"
    ac_mini=$(ac_json "$M/run"); require_json_values "$ac_mini"
    {
        printf '{\n'
        printf '  "frames": %s,\n' "$NMINI"
        printf '  "collisions": %s,\n' "$col_mini"
        printf '  "accommodation_coefficients_bottom_wall": %s\n' "$ac_mini"
        printf '}\n'
    } > "$OUT/mini/expected.json"

    # provenance (changes on every run: date, host load)
    {
        printf 'date: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        printf 'uname: %s\n' "$(uname -srm)"
        printf 'v2 checkout: %s\n' "$(git -C "$V2" describe --tags --always 2>/dev/null || echo unknown) $(git -C "$V2" rev-parse HEAD 2>/dev/null || true)"
        printf 'np: %s\nOMP_NUM_THREADS: %s\n' "$NP" "$OMP_NUM_THREADS"
        printf -- '--- %s -h | head -3\n' "$LMP"; "$LMP" -h | head -3
        printf -- '--- %s --version\n' "$MPIRUN"; "$MPIRUN" --version 2>&1 | head -1
        printf -- '--- %s --version\n' "$CXX"; "$CXX" --version 2>&1 | head -2
        printf -- '--- processor grid\n'; grep -h 'MPI processor grid' "$RUN/md/log.equil" "$RUN/md/log.meas" 2>/dev/null || true
        printf -- '--- md wall times\n'; cat "$RUN/md/md_times.txt" 2>/dev/null || true
        printf -- '--- tool build status (name rc warnings)\n'; cat "$BIN/build_status.tsv" 2>/dev/null || true
        printf -- '--- v2.0.0 sources used (sha256)\n'
        ( cd "$V2" && sha initialisation/TypeOfWalls/Explicit_layered/combine.cpp \
            initialisation/TypeOfWalls/Explicit_layered/unitCell.dat \
            initialisation/TypeOfWalls/Explicit_layered/Parameters.dat \
            LAMMPS_input/Explicit_layered/in.equil LAMMPS_input/Explicit_layered/in.meas tools/*.cpp )
    } > "$OUT/provenance.txt"
fi

log "done"
