# Reference answers

These files were made with the code at tag `v2.0.0`, before any refactoring.
Every later change is checked against them: results must be identical, or
within the tolerance stated in each case's README.

| Directory | What it covers | Made from |
|---|---|---|
| `R1_beam_mini/` | Molecular-beam line (Ar on Pt): the frozen `LAMMPS_input/Beam_ArPt/in.beam` on a 1600-atom cut of the slab, three incident conditions plus one regression-only case, analysed with the frozen `extract_alpha.py` | `make.sh` in that directory |
| `R2_examples_mini/` | One case from `examples/`: a short copy of `examples/Kerogen` (Fuel 2022), and the output of all 21 programs in `tools/` on its 11-column dump. The programs misread that dump (known issue KI-8), so this set checks that a rebuilt tool behaves exactly as before, not the physics | `make.sh` in that directory |
| `R3_channel_mini/` | Channel set-up of Chen (2024), Chapter 3 (Ar between two explicit Pt walls, 300 K, equilibrium MD): a small copy of `LAMMPS_input/Explicit_layered`, and the output of all 21 programs in `tools/` | `make.sh` in that directory |
| `FROZEN.sha256` | sha256 of every file used by the published papers (inputs, sources, tables) | `python3 tests/regress/check_frozen.py --write` |

Both `make.sh` scripts take the path of a `v2.0.0` checkout, for example

```bash
git worktree add --detach ../sk-v2 v2.0.0
bash tests/reference/R3_channel_mini/make.sh ../sk-v2
```

## Where the thesis numbers come from

`R3_channel_mini` follows the set-up of Chapter 3 of

Chen, Y., 2024. *Molecular dynamics modelling of gas–surface interactions*.
PhD thesis, The University of Edinburgh.
[doi:10.7488/era/5430](https://doi.org/10.7488/era/5430)

"Fig. 3.4", "eq. 3.6", "eq. 3.7" and "section 3.3.1" in these files refer to
that thesis. Fig. 3.4(a–c) gives TMAC 0.49, NMAC 0.67, NEAC 0.63 and EAC 0.26
for Ar on Pt at 300 K (least squares over all collisions).

## LAMMPS build

Both cases were run with LAMMPS `stable_22Jul2025_update4` (commit 611ca3b),
built with CMake from `cmake/presets/basic.cmake` (KSPACE, MANYBODY, MOLECULE,
RIGID) plus EXTRA-FIX and REAXFF, `-D BUILD_MPI=on -D BUILD_OMP=off
-D BUILD_SHARED_LIBS=off -D CMAKE_BUILD_TYPE=Release -D FFT=FFTW3`, with Apple
clang 21.0.0 (C++17) and Open MPI 5.0.9 on macOS arm64.

`lmp -h` reports the build as `-modified` because the source tree used for it
had a few extra fix styles and one edited ReaxFF file. None of the
`LAMMPS_src/` styles of this repository were compiled in, and neither case
uses any non-standard style: both need only `atom_style full` (MOLECULE),
`pair_style lj/cut` and standard fixes and computes. A stock build of the
same version with MOLECULE should therefore run both cases. Byte-identical
dumps need the same build, compiler, MPI and number of ranks; otherwise
compare with the tolerances in each case's README.

## Checks

```bash
python3 tests/regress/check_frozen.py      # frozen paper files unchanged, none added
python3 tests/regress/check_hygiene.py     # no large data or model files tracked
python3 tests/regress/check_reference_changes.py origin/main   # CHANGES.md updated with any reference change
(cd tests/reference/R3_channel_mini/mini && shasum -a 256 -c SHA256SUMS)
```

## Changing a reference

Reference files are only regenerated with their `make.sh`, and every change
gets an entry in `CHANGES.md` saying what changed and why. The same holds for
`FROZEN.sha256`.
