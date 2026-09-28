# Edits to the v2.0.0 files

The three files in `inputs/` are byte-for-byte copies of the v2.0.0 files
below, with only the lines listed here replaced. `make.sh` makes these edits
itself (it checks that each replaced line is exactly the old text first) and
writes the plain `diff` output to `inputs/EDITS.diff`. Nothing else changes,
including the missing final newline in `Parameters.dat`.

`combine.cpp` and `unitCell.dat` are used unchanged from
`initialisation/TypeOfWalls/Explicit_layered/`.

## inputs/Parameters.dat

Source: `initialisation/TypeOfWalls/Explicit_layered/Parameters.dat`

| line | v2.0.0 | here | why |
|---|---|---|---|
| 1 | `0    200          # xLo,  xHi` | `0    40           # xLo,  xHi` | small lateral box; combine.cpp rounds 40 up to 11 unit cells = 43.12 Å |
| 2 | `0    200          # yLo,  yHi ` | `0    40           # yLo,  yHi ` | same, y (the trailing space is kept) |

Result: 4150 atoms (278 Ar + 3872 Pt), box 43.12 × 43.12 Å, channel height 30 Å.

## inputs/in.equil

Source: `LAMMPS_input/Explicit_layered/in.equil`

| line | v2.0.0 | here | why |
|---|---|---|---|
| 16 | `processors      8 8 2` | `processors      * * *` | let LAMMPS pick the grid for 4 ranks (8 8 2 needs 128) |
| 106 | `run             100000 every 100000 "write_restart restart.*"` | `run             50000 every 50000 "write_restart restart.*"` | shorter equilibration (50 ps) |

## inputs/in.meas

Source: `LAMMPS_input/Explicit_layered/in.meas`

| line | v2.0.0 | here | why |
|---|---|---|---|
| 16 | `processors      8 8 2` | `processors      * * *` | as above |
| 51 | `pair_coeff   1 1 ${epsilonGas}   ${sigmaGas}   ${rCut}` | `#pair_coeff   1 1 ${epsilonGas}   ${sigmaGas}   ${rCut}` | switch gas–gas interaction off ... |
| 52 | `#pair_coeff   1 1 0.0000 ${sigmaGas} 0.1` | `pair_coeff   1 1 0.0000 ${sigmaGas} 0.1` | ... using the line already in the deck (thesis 3.2.1) |
| 106 | `run             6000000 every 6000000 "write_restart restart.*"` | `run             230000 every 230000 "write_restart restart.*"` | 230 ps production, 1151 dump frames |

Why gas–gas off: with it on, the Ar mean free path (about 4 Å at this
density) is shorter than the 12 Å virtual plane, so most "collisions" the
tools count are Ar–Ar collisions and all four accommodation coefficients come
out near 0.95. The thesis Chapter-3 virtual-plane method assumes gas–gas
interactions are off. The frozen v2.0.0 deck has them on (to be recorded as a
known issue, KI-15 in the v3 plan).
