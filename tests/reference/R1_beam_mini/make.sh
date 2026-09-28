#!/usr/bin/env bash
# Regenerate the R1 reference answers: the TAML Ar-Pt molecular beam with the
# frozen v2.0.0 in.beam (not changed by a single byte) on a 10x10-cell slab,
# serial LAMMPS, then the frozen extract_alpha.py on each run.
#
# usage: make.sh [V2_CHECKOUT] [RUN_ROOT]
#   V2_CHECKOUT  checkout of tag v2.0.0 (default: the repository holding this script)
#   RUN_ROOT     scratch directory for the runs, about 50 MB (default: ${TMPDIR:-/tmp}/sk_R1_beam_mini)
# environment: LMP (default lmp), PYTHON (default python3)
#
# Takes about 6 minutes on one core. The stored files in this directory are
# overwritten; README.md and inputs.sha256 are not touched.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V2="$(cd "${1:-$HERE/../../..}" && pwd)"
ROOT="${2:-${TMPDIR:-/tmp}/sk_R1_beam_mini}"
LMP="${LMP:-lmp}"
PY="${PYTHON:-python3}"

N_INSERT=500
N_LOOPS=1
NCELL=10
MAX_BYTES=2000000        # every stored file stays below 2 MB
WALL_LEFT_MAX=10         # eps1_th30: at most 2 % of N_INSERT still near the wall at the end
DECK="$V2/LAMMPS_input/Beam_ArPt"
SUM="$HERE/summarise.py"

# stored name, run directory name (extract_alpha.py reads eps and theta from
# it), V_INC [m/s] from the eps -> V_INC table in Beam_ArPt/README.md,
# THETA [deg], TAILSTEPS.
# TAILSTEPS 12000 left 12 of 500 eps1_th30 atoms (2.4 %) in the wall region,
# 16000 leaves 5 (1.0 %). eps0.4_th75 keeps 12000 so that no atom has come
# back yet and extract_alpha.py takes its N_return < 30 branch (alpha = 1);
# the same case with 16000 returns 268 atoms and is stored as well.
CASES=(
  "eps1_th30             eps1_th30   2198 30 16000"
  "eps1_th45             eps1_th45   2198 45 16000"
  "eps0.4_th75           eps0.4_th75 1390 75 12000"
  "eps0.4_th75_tail16000 eps0.4_th75 1390 75 16000"
)

sha256() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$@"; else shasum -a 256 "$@"; fi
}

# the deck's own relative paths must resolve: runs/<case>/../../ArPt_slab.data
make_tree() {
  local tree="$1"
  rm -rf "$tree"
  mkdir -p "$tree/runs"
  cp "$DECK/in.beam" "$tree/in.beam"
  cp "$ROOT/slab_mini.data" "$tree/ArPt_slab.data"
  cmp "$DECK/in.beam" "$tree/in.beam"
}

run_case() {
  local tree="$1" run="$2" vinc="$3" theta="$4" tsteps="$5"
  local rdir="$tree/runs/$run"
  local args=(-in ../../in.beam -var V_INC "$vinc" -var THETA "$theta" -var PHI 0
              -var T_WALL 300 -var T_GAS 300
              -var n_insert "$N_INSERT" -var n_loops "$N_LOOPS" -var TAILSTEPS "$tsteps")
  mkdir -p "$rdir"
  echo "cd LAMMPS_input/Beam_ArPt/runs/$run && lmp ${args[*]}" > "$rdir/command.txt"
  local t0=$SECONDS
  (cd "$rdir" && "$LMP" "${args[@]}" > stdout.txt 2>&1)
  echo "$tree $run $((SECONDS - t0)) s" | tee -a "$ROOT/timings.txt"
  grep -q "^Simulation completed" "$rdir/log.beam"
}

store_case() {
  local rdir="$1" dst="$2"
  rm -rf "$dst"
  mkdir -p "$dst"
  "$PY" "$DECK/extract_alpha.py" "$rdir" --save-events > "$rdir/extract_alpha_stdout.txt"
  { echo "# sha256 size_bytes frames"; "$PY" "$SUM" dumpinfo "$rdir/dump_gas.lammpstrj"; } > "$dst/dump_info.txt"
  grep "^Total particles" "$rdir/log.beam" > "$dst/final_counts.txt"
  cp "$rdir/command.txt" "$rdir/extract_alpha_stdout.txt" "$dst/"
  cp "$rdir/DataDir/particle_stats.txt" "$dst/particle_stats.txt"
  if [ "$(wc -c < "$rdir/events_out.csv")" -lt "$MAX_BYTES" ]; then
    cp "$rdir/events_out.csv" "$dst/events_out.csv"
  else
    echo "WARNING: $dst/events_out.csv would be 2 MB or more, not stored" >&2
  fi
  echo "$(basename "$dst"): $(awk '/^Total particles in wall region:/ {print $NF}' "$rdir/log.beam")" \
       "of $N_INSERT atoms still in the wall region at the end"
}

# the eps1_th30 dump, gzipped; cut to its first frames only if it does not fit
store_th30_dump() {
  local rdir="$1" dst="$HERE/eps1_th30"
  "$PY" "$SUM" gzip "$rdir/dump_gas.lammpstrj" "$dst/dump_gas.lammpstrj.gz" "$MAX_BYTES" \
    > "$dst/dump_gz_info.txt"
  cat "$dst/dump_gz_info.txt"
  # extract_alpha.py on exactly the stored (possibly cut) dump
  local chk="$ROOT/stored_dump/eps1_th30"
  rm -rf "$chk"
  mkdir -p "$chk"
  gunzip -c "$dst/dump_gas.lammpstrj.gz" > "$chk/dump_gas.lammpstrj"
  "$PY" "$DECK/extract_alpha.py" "$chk" --save-events > "$dst/extract_alpha_stdout_stored_dump.txt"
  cmp "$chk/events_out.csv" "$dst/events_out.csv"
}

check_th30_tail() {
  local n
  n=$(awk '/^Total particles in wall region:/ {print $NF}' "$HERE/eps1_th30/final_counts.txt")
  if [ "$n" -gt "$WALL_LEFT_MAX" ]; then
    echo "WARNING: TAILSTEPS too short for eps1_th30 ($n > $WALL_LEFT_MAX atoms left near the wall)" >&2
  fi
}

check_determinism() {
  local tree="$ROOT/rerun_eps1_th30/LAMMPS_input/Beam_ArPt"
  make_tree "$tree"
  run_case "$tree" eps1_th30 2198 30 16000
  local s1 s2 verdict=DIFFERENT
  s1=$("$PY" "$SUM" dumpinfo "$ROOT/eps1_th30/LAMMPS_input/Beam_ArPt/runs/eps1_th30/dump_gas.lammpstrj")
  s2=$("$PY" "$SUM" dumpinfo "$tree/runs/eps1_th30/dump_gas.lammpstrj")
  [ "$s1" = "$s2" ] && verdict=IDENTICAL
  printf '# eps1_th30 dump, serial, fresh tree each: sha256 size_bytes frames\nrun1 %s\nrun2 %s\n%s\n' \
    "$s1" "$s2" "$verdict" | tee "$HERE/determinism.txt"
  [ "$verdict" = IDENTICAL ]
}

main() {
  mkdir -p "$ROOT"
  : > "$ROOT/timings.txt"
  (cd "$V2" && sha256 -c "$HERE/inputs.sha256")
  {
    date -u '+%Y-%m-%dT%H:%M:%SZ'
    uname -a
    "$LMP" -h | sed -n '2p;/^Compiler:/p;/^MPI v/p'
    "$PY" --version
    echo "N_INSERT=$N_INSERT N_LOOPS=$N_LOOPS NCELL=$NCELL np=1 (no mpirun)"
  } | tee "$ROOT/provenance.txt"

  "$PY" "$HERE/crop_slab.py" "$DECK/ArPt_slab.data" "$ROOT/slab_mini.data" "$NCELL"
  gzip -c -n -9 "$ROOT/slab_mini.data" > "$HERE/slab_mini.data.gz"
  (cd "$ROOT" && sha256 slab_mini.data) > "$HERE/slab_mini.sha256"

  local c key run vinc theta tsteps tree pairs=()
  for c in "${CASES[@]}"; do
    read -r key run vinc theta tsteps <<< "$c"
    tree="$ROOT/$key/LAMMPS_input/Beam_ArPt"
    make_tree "$tree"
    run_case "$tree" "$run" "$vinc" "$theta" "$tsteps"
    store_case "$tree/runs/$run" "$HERE/$key"
    pairs+=("$key=$HERE/$key")
  done
  check_th30_tail
  store_th30_dump "$ROOT/eps1_th30/LAMMPS_input/Beam_ArPt/runs/eps1_th30"

  "$PY" "$SUM" expected "$HERE/expected.json" \
    "meta.deck=LAMMPS_input/Beam_ArPt/in.beam (v2.0.0, unchanged)" \
    "meta.analysis=LAMMPS_input/Beam_ArPt/extract_alpha.py (v2.0.0, unchanged)" \
    "meta.slab_sha256=$(cut -d' ' -f1 "$HERE/slab_mini.sha256")" \
    "meta.np=1" "${pairs[@]}"

  check_determinism

  local big
  big=$(find "$HERE" -type f -size +"$((MAX_BYTES / 1024))"k)
  if [ -n "$big" ]; then
    echo "ERROR: stored files of 2 MB or more: $big" >&2
    exit 1
  fi
}

main "$@"
