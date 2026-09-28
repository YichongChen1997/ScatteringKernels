# R3: thesis Chapter-3 channel case, mini size (reference answers)

Reference answers for the Chapter-3 workflow: Ar gas between two explicit,
layered Pt walls, walls and gas at 300 K, equilibrium MD, gas-wall collisions
counted at a virtual plane 12 Å above the bottom surface. This is also the
beginner tutorial case (`cases/ch3_channel/` in the v3 plan).

Everything here was generated from a clean checkout of tag **v2.0.0**
(commit e4ad3ce737b2b0c60e890071f17c87be6479079d) by `make.sh`. Nothing in
the repository was changed to make it.

## The case

- Geometry and walls: `initialisation/TypeOfWalls/Explicit_layered/combine.cpp`
  with `unitCell.dat` unchanged and `Parameters.dat` edited to a 40 Å lateral
  box, which combine.cpp rounds up to 11 fcc cells = 43.12 × 43.12 Å.
  Channel height 30 Å, bottom surface at z = 0. 8 Pt layers per wall:
  1 surface layer (types 5 bottom / 6 top, 242 atoms each), 5 bulk layers
  (types 2 / 3, 1210 each), both thermostatted, and 2 frozen layers
  (type 4, 968 atoms for both walls). 4150 atoms in total, 278 of them Ar
  (type 1).
- MD decks: `LAMMPS_input/Explicit_layered/in.equil` and `in.meas`, edited
  only as listed in [EDITS.md](EDITS.md): `processors * * *`, equilibration
  50 000 steps (gas and walls NVT 300 K), production 230 000 steps (gas NVE,
  walls NVT), gas-gas interaction switched off in `in.meas`. Time step 1 fs,
  gas dump every 200 steps, so 1151 frames (steps 50 000 to 280 000).
- Post-processing: all 21 `tools/*.cpp`, each compiled with plain
  `g++ -o name name.cpp` (no flags, the same recipe as the v2.0.0 Makefile),
  run with this `Specification.dat`:

```
1151                  # nTimeSteps (frames to read)
1.0  200  0          # deltaT[fs]  tSkip(dump every)  skipTimeStep(frames)
12   30              # rCut (virtual plane z, A)  H [A]
300  300             # Tg  Tw [K]
39.948  195.084      # mG  mW [g/mol]
```

## How it was run

```
V2=<clean checkout of v2.0.0>
RUN=<scratch dir for big files, ~2 GB incl. Contours output; default $TMPDIR/sk_R3_run>
tests/reference/R3_channel_mini/make.sh "$V2" "$RUN"
```

(The stored files were made with `V2` = a clean detached checkout of tag
v2.0.0 and `RUN` in a scratch directory outside the repository.)

`make.sh` runs, in order (set `STAGES` to run only some of them):

1. `init`: copy combine.cpp and unitCell.dat, make Parameters.dat, then
   `g++ -o combine combine.cpp && ./combine` -> `data.dat`.
2. `md`: make in.equil / in.meas, then
   `mpirun -np 4 lmp -in in.equil` and `mpirun -np 4 lmp -in in.meas`
   (OMP_NUM_THREADS=1; LAMMPS picked a 2 by 1 by 2 processor grid).
3. `bin`: `g++ -o <name> <V2>/tools/<name>.cpp` for all 21 tools.
4. `full`: every tool on the full dump, each in its own directory with links
   to `dump_meas_gas.lammpstrj`, `dump_meas_wall.lammpstrj` and
   `Specification.dat`; stdin from /dev/null; exit code, wall time and peak
   memory recorded; a tool is killed above 8 GB resident memory or 3600 s.
5. `mini`: the first 80 frames of the gas dump, `gzip -9n`, then gunzip it
   again and run every tool on it with `nTimeSteps = 80`.
6. `collect`: copy the small files below into this directory.

It needs bash, perl, gzip, awk, pgrep/ps, and `sha256sum` or `shasum`.
`/usr/bin/time -l` (macOS) or `-v` (GNU) is used for peak memory if present.
With the defaults, `collect` overwrites the stored files in place, so after a
rerun `git diff` shows what changed. Only `provenance.txt` and the `wall_s`
and `max_rss_MB` columns of `tool_status.tsv` should change.

## Environment (2026-09-28)

- macOS (Darwin 25.5.0, arm64), 10 cores, 32 GB, shared with other jobs
- `lmp -h`: `Large-scale Atomic/Molecular Massively Parallel Simulator - 22 Jul 2025 - Update 4`,
  `Git info (stable / stable_22Jul2025_update4-modified)`; how it was built is
  in `../README.md` (LAMMPS build)
- `mpirun (Open MPI) 5.0.9`, **np = 4** (use the same np to get the same bytes)
- `g++` = Apple clang version 21.0.0 (clang-2100.1.1.101), target arm64-apple-darwin25.5.0,
  no flags (clang's default language standard)
- Full details, including the sha256 of every v2.0.0 source file used, are
  in `provenance.txt`.

## Wall-clock times (this machine, under load)

| step | time |
|---|---|
| combine | < 1 s |
| in.equil, 50 000 steps, 4 ranks | 63 s |
| in.meas, 230 000 steps, 4 ranks | 289 s (an earlier identical run under heavier load: 507 s) |
| compiling 21 tools | about 10 s |
| each tool on the full dump | 1.3-2.0 s, except pp_meas_Contours 46 s and pp_meas_wallTemp 0.5 s |
| each tool on the 80-frame dump | 0.26 s, except pp_meas_Contours 44 s |
| whole make.sh, all stages, from scratch | 465 s (equil 58 s, meas 268 s in that run) |

## Results on the full run (`full/expected.json`)

Bottom wall, all collisions (pp_meas_ACs; top wall is not counted by the
tool). 1981 collisions started, 1869 ended and were used.

| | general formula (eq. 3.6) | least squares (eq. 3.7) | thesis Fig. 3.4 |
|---|---|---|---|
| TMAC | 0.535549 | 0.495551 | 0.49 |
| NMAC | 0.378085 | 0.670302 | 0.67 |
| NEAC | 0.345013 | 0.648138 | 0.63 |
| EAC | 0.156304 | 0.264906 | 0.26 |
| (AlphaEx, tangential-x energy) | 0.289702 | 0.569353 | - |

Where these come from: in every `Alpha*.txt`, column 1 is the velocity
window, columns 2-3 are the windowed ("partial") general / least-squares
values, and columns 4-5 are the values over all collisions, general formula
then least squares. Columns 4-5 are the same on every row; `make.sh` checks
this. The least-squares column matches Fig. 3.4. The general formula is
poor here because at equilibrium its denominator is close to zero (the
limitation discussed in thesis 3.3.1).

sha256 of the big files (not stored):

- `data.dat` (125 259 bytes): `5beb62ff3ceec51742597aa0867a4d67699bac0766e2555d0b10327772da8b9f`
- `dump_meas_gas.lammpstrj` (43 718 663 bytes, 1151 frames, 16 columns): `671231ee51353d8bd5125f92a704a14eb2e64d51ff80cb5770d5bfcc951fd69b`
- `dump_meas_wall.lammpstrj` (195 155 bytes, 1 frame at step 200 000): `b01b2b43a77e585189feb537622aadf0bf97ed4f3059f84e35385e07cebf148c`

`full/OUTPUT_SHA256SUMS` has the sha256 of every tool output on the full
dump (the outputs themselves are not stored).

## The byte-exact regression set (`mini/`)

- `dump_meas_gas.first80.lammpstrj.gz`: the first 80 frames (steps 50 000
  to 65 800) of the gas dump, `gzip -9n` (1 330 413 bytes).
- `dump_meas_wall.lammpstrj.gz`: the whole wall dump (74 748 bytes).
- `Specification.dat`: as above with `nTimeSteps = 80`.
- `outputs/<tool>/`: every file each tool wrote, plus its stdout
  (`_stdout.txt`) and non-empty stderr (`_stderr.txt`). Files of 2 MB or
  more are not stored; `outputs/NOT_STORED.txt` gives their sha256 and size
  (only pp_meas_Contours `PES.txt` and `Fz.txt`, about 307 MB each).
- `SHA256SUMS`: the three inputs and every stored output, so that
  `(cd mini && shasum -a 256 -c SHA256SUMS)` passes as it is.
- `DUMP_CONTENT.sha256`: sha256 of the two dumps after `gunzip`. The `.gz`
  bytes depend on the gzip build, so compare regenerated dumps by content.
- `expected.json`: collisions and ACs for the 80 frames. Only 54 collisions
  end in 80 frames, so these numbers are for byte comparison only and have no
  physical meaning (NMAC general formula is even negative).

To check a build: gunzip both dumps into an empty directory, add
`Specification.dat`, run each tool there, and compare its outputs with the
`outputs/` lines of `SHA256SUMS` (and `NOT_STORED.txt` for the two large
files).

## Tool status

Same on the full and the 80-frame set (`full/tool_status.tsv`,
`mini/tool_status.tsv`, with stdout tails in `tool_stdout_tail.txt`):

- 19 tools exit 0 with normal output.
- **pp_meas_Traj crashes** (exit 139, segmentation fault). It tracks
  hard-coded atom ids 53, 235, 1000 and 5000; this case has only 278 gas
  atoms, so reading the trajectory of id 1000 goes past the end of an empty
  vector. It writes a partial `molecules.trj` first. Recorded, not fixed. The
  partial output was the same in every run here.
- **pp_meas_Contours** needs 4.2 GB peak memory (800 × 800 × 240 bins,
  four arrays) and writes two text files of about 307 MB, whatever the
  number of frames. It stayed under the 8 GB limit and finished.
- **pp_meas_wallTemp** exits 0 but gives no useful numbers. It reads
  `nTimeSteps / 100` frames from the wall dump and skips frame 0 when
  writing. The wall dump here has one frame (step 200 000), so on the full
  run it writes 10 rows of zeros (frames 1-10 are read past the end of the
  file) and on the 80-frame set it reads 0 frames and writes empty files.
- Other outputs the v3 plan already marks as buggy (pp_meas_VDFF `xi_Mag*`,
  pp_meas_VDFPi `xi_TxPrime_beam.txt`) are stored as they are; they are
  regression data, not physics.

## Things that look odd but are expected

- **Thermo shows T ≈ 230 K, not 300 K.** LAMMPS `thermo temp` averages over
  all 4150 atoms, including the 968 frozen Pt atoms that never move.
  300 K × (4150 − 968) / 4150 = 230 K. The mean thermo T was 230.0 K in
  equilibration (after the first 10 ps) and 231.3 K in production.
- `in.meas` says `boundary p p p`, but the dump header says `pp pp ff`:
  `read_restart` takes the boundary from the restart file (`p p f` from
  in.equil).
- The wall dump has 1 frame, not 2: production starts at step 50 000 and
  the wall dump interval is 200 000, and LAMMPS writes only at multiples of
  the interval.
- pp_meas_ACs prints "Collisions Started" (1981) > "Collisions Ended" (1869):
  a collision counts only once the molecule leaves the 12 Å layer again.
  The 112 missing are exactly the molecules still inside z < 12 Å in the
  last frame. Molecules already inside in the first frame are not counted
  as started.
- The frozen v2.0.0 in.meas has gas-gas interaction on. With it on, all four
  ACs come out near 0.95 (the mean free path, about 4 Å, is shorter than the
  12 Å virtual plane, so Ar-Ar collisions are counted). See EDITS.md.

## Tolerance on another machine

A different LAMMPS build, compiler, MPI or number of ranks gives different
dump bytes, so the full-run numbers are then compared with a tolerance, not
byte for byte. The statistical error of each least-squares coefficient was
estimated by bootstrap (400 resamples of the 1869 bottom-wall collisions
written by `pp_meas_Correlation`):

| coefficient | stored value | bootstrap sd | tolerance, 3 x sqrt(2) x sd |
|---|---|---|---|
| TMAC | 0.495551 | 0.020 | 0.09 |
| NMAC | 0.670302 | 0.022 | 0.09 |
| NEAC | 0.648138 | 0.027 | 0.11 |
| EAC | 0.264906 | 0.021 | 0.09 |

The sqrt(2) is there because both runs carry their own noise. The general
formula (column 4) is not used for comparisons: in equilibrium its
denominator is close to zero and it is not a stable number.

## Is a rerun byte-identical?

On this machine, yes, with the same LAMMPS binary and np = 4:

- A second run of the whole `make.sh`, from scratch into a new directory,
  reproduced every stored file byte for byte (only `provenance.txt` and the
  `wall_s` and `max_rss_MB` columns of `tool_status.tsv` differed).
- This run and an earlier, separate run of the same inputs (about 45 min
  apart, different machine load) gave byte-identical `restart.50000`,
  `dump_equil_gas.lammpstrj`, gas dump and wall dump.
- The tools gave byte-identical outputs on repeated runs, including the
  partial output of the crashing pp_meas_Traj.
- `combine` calls `srand(time(NULL))` but never calls `rand()` (its
  `fRand` is unused), and the missing `roughSurface.dat` just gives a flat
  surface, so `data.dat` is the same every time.
- `gzip -n` leaves out the file name and time stamp, so the `.gz` files are
  reproducible (with the same gzip).

Expect different dump bytes with another np, another LAMMPS build or version,
another MPI, or another CPU (floating-point summation order changes and the
trajectories drift apart). Then compare ACs with a tolerance, not bytes. The
tool outputs for a given dump should be the same on any machine only if the
compiler and C++ library print floating-point numbers the same way; check
that with the `mini/` set.

## Files

```
README.md            this file
EDITS.md             every line changed from v2.0.0, and why
make.sh              regenerates everything
provenance.txt       date, versions, np, sha256 of the v2.0.0 sources used
inputs/              Parameters.dat, in.equil, in.meas exactly as run; EDITS.diff
full/                Specification.dat, SHA256SUMS (data.dat, dumps),
                     OUTPUT_SHA256SUMS (tool outputs), expected.json,
                     tool_status.tsv, tool_stdout_tail.txt
mini/                the 80-frame regression set (see above)
```

Note: the repository `.gitignore` ignores `*.txt`, `*.out`, `data.dat`,
`*.lammpstrj` and `restart.*`; several files here match `*.txt`, so an
exception for `tests/reference/**` is needed before committing.
