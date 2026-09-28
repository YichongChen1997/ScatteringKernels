# R1: Ar-Pt molecular beam on a small slab (TAML line)

Reference answers for the molecular-beam line, made with the v2.0.0 deck
`LAMMPS_input/Beam_ArPt/in.beam` unchanged (sha256 in `inputs.sha256`) and the
v2.0.0 `extract_alpha.py` unchanged. The only thing that differs from a
production run is the slab: `ArPt_slab.data` is cropped from 103 x 103 to
10 x 10 fcc cells (39.2 x 39.2 Angstrom, all 8 layers, 1600 Pt, 400 of them
type 3 fixed), and 500 atoms are fired instead of 2500.

## Regenerate

    tests/reference/R1_beam_mini/make.sh <v2.0.0 checkout> <scratch dir>

Used here:

    make.sh ../sk-v2 /tmp/sk-ref-R1

`make.sh` checks the three frozen inputs against `inputs.sha256`, crops the
slab with `crop_slab.py`, builds one tree per case so that the deck's own
relative path resolves,

    <scratch>/<case>/LAMMPS_input/Beam_ArPt/in.beam            copy of the frozen deck
    <scratch>/<case>/LAMMPS_input/Beam_ArPt/ArPt_slab.data     the 1600-atom slab
    <scratch>/<case>/LAMMPS_input/Beam_ArPt/runs/eps<E>_th<T>/ run directory

runs LAMMPS in serial, runs `extract_alpha.py <rundir> --save-events`, copies the
small outputs here, runs eps1_th30 a second time in a fresh tree and compares
the dump sha256. Nothing is written into the V2 checkout. About 6 minutes.

## Environment of the stored run

- Date: 2026-09-28, 13:19 UTC.
- Machine: Apple T6000 (arm64), Darwin 25.5.0, heavily loaded by other jobs.
- LAMMPS: `Large-scale Atomic/Molecular Massively Parallel Simulator - 22 Jul 2025 - Update 4`
  (git `stable_22Jul2025_update4-modified`); how it was built is in
  `../README.md` (LAMMPS build).
- Compiler: Clang C++ Apple LLVM 21.0.0 (clang-2100.1.1.101), OpenMP not enabled, C++17.
- MPI library: Open MPI 5.0.9, but np = 1: `lmp` is started without `mpirun`.
- Python 3.13.3 for `crop_slab.py`, `extract_alpha.py` and `summarise.py`.

## Cases

Each run directory gets exactly this command (also in `<case>/command.txt`):

    cd LAMMPS_input/Beam_ArPt/runs/eps1_th30   && lmp -in ../../in.beam -var V_INC 2198 -var THETA 30 -var PHI 0 -var T_WALL 300 -var T_GAS 300 -var n_insert 500 -var n_loops 1 -var TAILSTEPS 16000
    cd LAMMPS_input/Beam_ArPt/runs/eps1_th45   && lmp -in ../../in.beam -var V_INC 2198 -var THETA 45 -var PHI 0 -var T_WALL 300 -var T_GAS 300 -var n_insert 500 -var n_loops 1 -var TAILSTEPS 16000
    cd LAMMPS_input/Beam_ArPt/runs/eps0.4_th75 && lmp -in ../../in.beam -var V_INC 1390 -var THETA 75 -var PHI 0 -var T_WALL 300 -var T_GAS 300 -var n_insert 500 -var n_loops 1 -var TAILSTEPS 12000
    cd LAMMPS_input/Beam_ArPt/runs/eps0.4_th75 && lmp -in ../../in.beam -var V_INC 1390 -var THETA 75 -var PHI 0 -var T_WALL 300 -var T_GAS 300 -var n_insert 500 -var n_loops 1 -var TAILSTEPS 16000

V_INC comes from the eps table in `Beam_ArPt/README.md`. The wall
equilibration is the deck's fixed 10000 steps; the seeds are the deck defaults.

| stored as | TAILSTEPS | frames | left near wall | events | alpha_t | alpha_n | eac | LAMMPS wall time |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `eps1_th30` | 16000 | 525 | 5 | 495 | 0.4646 +- 0.0266 | 0.7392 +- 0.0096 | 0.5915 | 73 s |
| `eps1_th45` | 16000 | 525 | 7 | 493 | 0.3343 +- 0.0166 | 0.5597 +- 0.0167 | 0.4660 | 74 s |
| `eps0.4_th75` | 12000 | 445 | 500 | 0 | 1 (trap-limit branch) | 1 | 1 | 59 s |
| `eps0.4_th75_tail16000` | 16000 | 525 | 232 | 268 | 0.2189 +- 0.0076 | -111.67 +- 5.16 | 0.0807 | 74 s |

Full-precision values are in `expected.json`. "Left near wall" is the
`Total particles in wall region` line of `log.beam` (0 < z < 15 Angstrom).

**`eps0.4_th75_tail16000` is for regression only.** At the end of its tail 232 of
the 500 atoms (46%) are still near the wall, so its coefficients come from an
incomplete tail and are not a physical result. It is kept because it is the
only case that sends 75-degree atoms through the measured branch of
`extract_alpha.py`.

## Tolerance on another machine

Serial runs with the same `lmp` binary are byte-identical (see
`determinism.txt`). With a different LAMMPS build, compiler or number of
ranks the dumps differ, and the coefficients are compared with
|delta alpha| <= 3 x sqrt(err_a^2 + err_b^2), where err are the bootstrap
errors printed by `extract_alpha.py` for the two runs. The event count may
differ by a few atoms. Compare `.gz` files by the sha256 of their uncompressed
content, since gzip output depends on the zlib build.

## Stored files

- `make.sh`, `crop_slab.py`, `summarise.py`: the generator.
- `inputs.sha256`: sha256 of the frozen `in.beam`, `extract_alpha.py`,
  `ArPt_slab.data` (paths relative to the V2 checkout).
- `slab_mini.data.gz`, `slab_mini.sha256`: the cropped slab and the sha256 of
  its uncompressed text.
- `expected.json`: alpha_t, alpha_n, eac, their bootstrap errors, events, and
  the run variables, per case.
- `determinism.txt`: sha256, size and frame count of the eps1_th30 dump from
  two independent runs.
- Per case: `command.txt`, `dump_info.txt` (sha256, size, frame count of
  `dump_gas.lammpstrj`), `extract_alpha_stdout.txt`, `events_out.csv`,
  `final_counts.txt` (the `Total particles` lines of `log.beam`),
  `particle_stats.txt` (the deck's `DataDir/particle_stats.txt`).
- `eps1_th30/dump_gas.lammpstrj.gz`: the whole eps1_th30 dump (1.40 MB gzip,
  5.53 MB text, all 525 frames, not cut). `dump_gz_info.txt` gives the sha256
  of the uncompressed text, which equals the one in `dump_info.txt`, and
  `extract_alpha_stdout_stored_dump.txt` is `extract_alpha.py` run on the
  unpacked file.

The repository `.gitignore` ignores `*.txt`; these files need the
`!tests/reference/**` exception.

## Notes

- **TAILSTEPS.** 12000 left 12 of 500 eps1_th30 atoms (2.4 %) and 16 of 500
  eps1_th45 atoms in the wall region, above the 2 % limit, so both use 16000
  (5 and 7 left). With 12000 the results were alpha_t = 0.4752, alpha_n = 0.7351,
  events 488 (th30) and 0.3428, 0.5514, events 484 (th45), inside one error bar
  of the 16000 values.
- **eps0.4_th75 reaches the alpha = 1 branch only because the run stops too early.**
  The normal speed is 1390 cos 75 = 360 m/s, so the atoms need about 10 ps to
  fall from the insertion plane (z = 49-50) to the force cutoff at z = 15. At the
  end of the 12000-step tail all 500 are still falling (z = 5.0 to 6.0, every
  vz < 0) and none has touched the surface, so `extract_alpha.py` finds 0 returns
  and sets alpha_t = alpha_n = eac = 1. This case is kept at 12000 on purpose, as
  the regression case for that branch. With 16000, 268 atoms come back and the
  measured branch runs; it is stored as `eps0.4_th75_tail16000`. The same slow
  approach may be behind the "no atom returns at 75 degrees" statement in
  `Beam_ArPt/README.md` (default TAILSTEPS 8000); this was not checked on the full slab.
- **alpha_n = -111.7 at eps0.4_th75_tail16000** is the frozen formula, not a run
  problem: the incident normal energy is 1.036 kT_w, so the denominator
  (E_n,i - 1) is 0.036.
- The `[trap-limit]` line of `extract_alpha.py` prints `1 - N_return/2500` with
  2500 fixed in the code; with 500 inserted atoms that fraction is only right
  because N_return = 0.
- **Determinism.** Two serial runs of eps1_th30 in fresh trees gave the same dump
  sha256 (`80ceda2a...0633`). The whole set was generated three times (tails of
  12000 and 16000, then the final `make.sh`), and every overlapping stored file
  matched byte for byte. The gzip files are byte-identical for the same gzip
  (`slab_mini.data.gz`, macOS `gzip -n -9`) or zlib (the dump, Python `gzip`,
  mtime 0) build; compare the uncompressed sha256 across machines. A different
  LAMMPS build, compiler or process count is not expected to reproduce the dump
  bytes.
