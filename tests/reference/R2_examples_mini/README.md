# R2: examples/Kerogen, mini size (reference answers)

Reference answers for one case in `examples/`: methane flowing between two
rough kerogen surfaces, driven by a constant force, at 423 K
(`examples/Kerogen`). This is the example that goes with

Chen, Y., Li, J., Datta, S., Docherty, S.Y., Gibelli, L. and Borg, M.K., 2022.
Methane scattering on porous kerogen surfaces and its impact on mesopore
transport in shale. *Fuel*, 316, p.123259.

It is a shortened test copy, not a reproduction of the paper.

**Why this example.** Of the three cases in `examples/`, Kerogen is the only
one that runs in minutes on a laptop: its `data.dat` has 6189 atoms, against
about 3 MB for `Contaminated_adsorbate` and 54 MB for `Fixed_adsorbate`
(whose decks also ask for 512 MPI ranks). Its `combine.cpp` is used
unchanged, so the structure is `EFK_50A_0.80.xyz`, the file it opens by name.

**What R2 is for.** The Kerogen decks write an 11-column gas dump
(`id type x y z vx vy vz` and three stress components). All 21 programs in
`tools/` read 16 values per atom, so they misread this dump without any
error message (known issue KI-8). R2 stores exactly what each program
writes when given this real `examples/` dump. It is a byte-for-byte check
that a rebuilt tool behaves as the v2.0.0 one did, not a physical result.
The channel case `../R3_channel_mini/` (16-column dump) is the reference
that exercises the tools' physics.

Everything here was generated from a clean checkout of tag **v2.0.0**
(commit e4ad3ce737b2b0c60e890071f17c87be6479079d) by `make.sh`. Nothing in
the repository was changed to make it.

## The case

- Initial configuration: `examples/Kerogen/initialisation/combine.cpp`,
  unchanged. It opens `EFK_50A_0.80.xyz` by name (the only file it reads)
  and writes `data.dat`: box 50 × 50 Å in x and y, z from −32 to 132 Å,
  channel height H = 100 Å, 6189 atoms: 87 methane molecules (type 1, one
  LJ site each, 2 MPa), two fixed barrier walls of 625 atoms (type 2), and
  two rough kerogen slabs of 2426 atoms each (C, O, H = types 3, 4, 5).
  `combine.cpp` calls `srand` but never `rand`, so `data.dat` is the same
  on every run.
- MD decks: `examples/Kerogen/in.equil` and `in.meas`, edited only as
  listed in [EDITS.md](EDITS.md): `processors * * *`, equilibration
  20 000 steps (10 ps), production 100 000 steps (50 ps). Time step 0.5 fs,
  gas dump every 500 steps, so 201 frames (steps 20 000 to 120 000).
- Post-processing: all 21 `tools/*.cpp`, each compiled with plain
  `g++ -o name name.cpp` (no flags, as in the v2.0.0 Makefile), run on the
  whole production dump with this hand-written `Specification.dat`:

```
201                  # nTimeSteps (frames to read)
0.5  500  0          # deltaT[fs]  tSkip(dump every)  skipTimeStep(frames)
15   100             # rCut (virtual plane z, A)  H [A]
423  423             # Tg  Tw [K]
16.043  12.0107      # mG  mW [g/mol]
```

  The virtual plane is put at the 15 Å gas–wall cut-off of the decks; mW
  (kerogen carbon) is only used by `pp_meas_wallTemp`. The whole dump fits
  in 0.67 MB gzip'd, so nothing was cut: the regression dump
  `mini/dump_meas_gas.first201.lammpstrj.gz` is the full production dump.

## How it was run

```
V2=<clean checkout of v2.0.0>
RUN=<scratch dir for big files, about 1.6 GB incl. Contours output; default $TMPDIR/sk_R2_run>
LMP=<path to lmp> tests/reference/R2_examples_mini/make.sh "$V2" "$RUN"
```

`make.sh` runs, in order (set `STAGES` to run only some of them):

1. `init`: copy combine.cpp and EFK_50A_0.80.xyz, `g++ -o combine combine.cpp && ./combine` -> `data.dat`.
2. `md`: make in.equil / in.meas, then `mpirun -np 4 lmp -in in.equil` and
   `mpirun -np 4 lmp -in in.meas` (OMP_NUM_THREADS=1; LAMMPS picked a
   1 by 1 by 4 processor grid).
3. `bin`: `g++ -o <name> <V2>/tools/<name>.cpp` for all 21 tools (no warnings).
4. `tools`: `gzip -9n` the first 201 frames of the gas dump, gunzip it
   again and run every tool on it, each in its own directory with links to
   `dump_meas_gas.lammpstrj` and `Specification.dat`; stdin from /dev/null;
   exit code, wall time and peak memory recorded; a tool is killed above
   8 GB resident memory or `TOOL_TIMEOUT` (600 s).
5. `repeat`: the same again in a fresh directory; stops with an error
   unless every output is byte-identical to stage 4.
6. `collect`: copy the small files below into this directory.

It needs bash (3.2 is enough), perl, gzip, awk, pgrep/ps, and `sha256sum`
or `shasum`. `/usr/bin/time -l` (macOS) or `-v` (GNU) is used for peak
memory if present. With the defaults, `collect` overwrites the stored
files in place, so after a rerun `git diff` shows what changed. Only
`provenance.txt` and the `wall_s` and `max_rss_MB` columns of
`mini/tool_status.tsv` should change.

## Environment (2026-09-30)

- macOS (Darwin 25.5.0, arm64), 10 cores, 32 GB, shared with other jobs
  (load average 8 to 14 during the run)
- `lmp -h`: `Large-scale Atomic/Molecular Massively Parallel Simulator - 22 Jul 2025 - Update 4`,
  `Git info (stable / stable_22Jul2025_update4-modified)`; how it was built
  is in `../README.md` (LAMMPS build). The decks need only MOLECULE
  (`atom_style full`), `pair_style lj/cut` and standard fixes and computes.
- `mpirun (Open MPI) 5.0.9`, **np = 4** (use the same np to get the same bytes)
- `g++` = Apple clang version 21.0.0 (clang-2100.1.1.101), target arm64-apple-darwin25.5.0,
  no flags
- Full details, including the sha256 of every v2.0.0 source file used, are
  in `provenance.txt`.

## Wall-clock times (this machine, under load)

| step | time |
|---|---|
| combine | < 1 s |
| in.equil, 20 000 steps, 4 ranks | 57 s |
| in.meas, 100 000 steps, 4 ranks | 289 s |
| compiling 21 tools | about 9 s |
| each tool on the 201-frame dump | 0.5-1.0 s, except pp_meas_Contours 55 s (5.3 GB peak memory) |
| whole make.sh, all stages, from scratch | 485 s (equil 56 s, meas 285 s in that run) |

## What the tools did (KI-8)

Every tool reads the frame header correctly, then reads 16 values per atom
with `>>`. For the first molecule of the first frame it takes the 11 values
of that line plus the first five of the next line (id, type, x, y, z of
molecule 2) as KE, PE, fx, fy, fz. The next read wants an integer id and
finds `-0.000509214` (vx of molecule 2): it takes `-0`, and the read of the
type then fails on the `.`. From there the input stream is in a failed
state, every later read does nothing, and the tools run through all 201
frames without an error, printing `currentTimeStep = 20000` for each one.
The failed read sets the type to 0; every read after it leaves its
variable as it was, so every later atom slot holds molecule 1's values.
In effect the tools see molecule 1, at z = 8.5 Å, and nothing else:

- The tools that follow gas molecules by type (ACs, Collisions,
  Correlation, penACs, penBeam, penOverall, AngularF, AngularP, FAngular,
  PAngular, ResidenceTime, Displacement, VDFF, VDFPi, VDFPr) find no
  collision: the counts are 0 and every output column holds 0 or `nan`.
  `pp_meas_Correlation/bottomVelocities.txt` is empty.
- `pp_meas_Bins` and `pp_meas_GasOverTime` do not check the type, so they
  count all 87 slots as copies of molecule 1 in every frame: Bins has one
  non-empty bin, at z = 8.5 Å, and gasOverTime.txt has 201 identical rows
  (mean velocity 855.857 m/s is molecule 1's vx).
- `pp_meas_Traj` crashes with signal 11 (exit 139) after writing only the
  9-line header of `molecules.trj` (KI-7), as in R3.
- `pp_meas_NearWall` writes a 12-byte `nearWallMolecules.xyz` (a count of
  0 and one line of zeros).
- `pp_meas_wallTemp`: the Kerogen decks write no wall dump, so the tool
  opens a file that does not exist and never sets `nAtoms` or
  `currentTimeStep`. It reads 201 / 100 = 2 "frames" (KI-6) and writes one
  row `2.5e-13 0 0 0` to each output file.
- `pp_meas_Contours` allocates 800 × 800 × 300 bins (300 = 20 × rCut),
  needs 5.3 GB of memory and writes two 384 MB text files whatever it
  reads. They are not stored; `mini/outputs/NOT_STORED.txt` has their
  sha256 and size.

Two outputs depend on variables the v2.0.0 code never initialises:
`pp_meas_wallTemp` (all of it) and the pressure column of
`pp_meas_GasOverTime/gasOverTime.txt` (`stressTemp`, KI-4). They were the
same in every run here, but a different compiler or build may change them.

## Is a rerun byte-identical?

On this machine, yes, with the same LAMMPS binary and np = 4:

- The `repeat` stage ran all 21 tools a second time in a fresh directory:
  all 66 files (every output and stdout, including the two Contours files
  and the partial output of the crashing pp_meas_Traj) were byte-identical.
- A second run of the whole `make.sh`, from scratch into a new run
  directory and a new output directory, reproduced every stored file byte
  for byte, and gave the same `data.dat`, gas dumps (`RUN_SHA256SUMS`) and
  restart files. Only `provenance.txt` and the `wall_s` and `max_rss_MB`
  columns of `mini/tool_status.tsv` differed.
- `gzip -n` leaves out the file name and time stamp, so the `.gz` file is
  reproducible with the same gzip. `mini/DUMP_CONTENT.sha256` has the
  sha256 of the uncompressed dump, which is what matters.

## Tolerance on another machine

No tolerance: this set is compared byte for byte. The tool outputs depend
only on the stored dump and `Specification.dat`, not on the MD. To check a
new build of the tools, run each one on the unpacked `mini/` dump (as
stage `tools` does) and compare every file with `mini/outputs/`, and the
two Contours files with the sha256 in `mini/outputs/NOT_STORED.txt`. With
another compiler or C++ library, a difference in how floating-point
numbers are printed, or in the two uninitialised outputs above, is
possible; anything else that differs is a real change in behaviour.

A new MD run with another LAMMPS build, MPI, CPU or number of ranks gives
different dump bytes (`RUN_SHA256SUMS` will not match). There is no
physical number to compare with a tolerance here; such a run is only
useful to see that the edited decks still run and give 201 frames of 87
molecules.

## Things that look odd but are expected

- **Thermo shows T ≈ 338 K, not 423 K.** LAMMPS `thermo temp` averages over
  all 6189 atoms, including the 1250 barrier atoms that never move:
  423 K × (6189 − 1250) / 6189 = 337.6 K. The mean over production was
  337.6 K.
- **PotEng is about 1.39e6 kcal/mol and Press about 1.15e6 atm.** The fixed
  barrier atoms sit 2 Å apart with the wall–wall LJ σ = 2.471 Å; they never
  move, so this is a constant offset.
- `in.meas` keeps methane–methane interaction on and the driving force;
  this is the v2.0.0 example as written.

## Files

```
README.md            this file
EDITS.md             every line changed from v2.0.0, and why
make.sh              regenerates everything
provenance.txt       date, versions, np, sha256 of the v2.0.0 sources used
RUN_SHA256SUMS       sha256 of init/data.dat and of the full equilibration
                     and production gas dumps (paths relative to RUN)
inputs/              in.equil, in.meas exactly as run; EDITS.diff
mini/                Specification.dat, the gzip'd 201-frame gas dump,
                     DUMP_CONTENT.sha256, SHA256SUMS (every stored file),
                     tool_status.tsv, tool_stdout_tail.txt, expected.json,
                     outputs/<tool>/ (every output < 2 MB, and the tool's
                     stdout), outputs/NOT_STORED.txt
```
