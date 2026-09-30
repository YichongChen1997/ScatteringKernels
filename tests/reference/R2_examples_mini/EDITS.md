# Edits to the v2.0.0 files

The two files in `inputs/` are byte-for-byte copies of
`examples/Kerogen/in.equil` and `examples/Kerogen/in.meas` at tag v2.0.0,
with only the lines listed here replaced. `make.sh` makes these edits itself
(it checks that each replaced line is exactly the old text first) and writes
the plain `diff` output to `inputs/EDITS.diff`. Nothing else changes.

`examples/Kerogen/initialisation/combine.cpp` and the structure it opens by
name, `examples/Kerogen/initialisation/EFK_50A_0.80.xyz`, are used unchanged.
No other file of the example is read.

## inputs/in.equil

Source: `examples/Kerogen/in.equil`

| line | v2.0.0 | here | why |
|---|---|---|---|
| 17 | `processors      16 16 2` | `processors      * * *` | let LAMMPS pick the grid for 4 ranks (16 16 2 needs 512) |
| 74 | `run             200000 every 100000 "write_restart restart.*"` | `run             20000 every 10000 "write_restart restart.*"` | 10 ps instead of 100 ps; `every` scaled by the same factor, so the run is still split in 2 parts |

## inputs/in.meas

Source: `examples/Kerogen/in.meas`

| line | v2.0.0 | here | why |
|---|---|---|---|
| 17 | `processors      16 16 2 ` | `processors      * * * ` | as above (the trailing space of the original line is kept) |
| 90 | `run             10000000 every 1000000 "write_restart restart.*"` | `run             100000 every 10000 "write_restart restart.*"` | 50 ps instead of 5 ns; `every` scaled by the same factor, so the run is still split in 10 parts |

With a time step of 0.5 fs and the gas dump every 500 steps, production
(steps 20 000 to 120 000) gives 201 dump frames.

Everything else is as in v2.0.0, including the settings that make this a
flow case rather than an equilibrium one: methane–methane interaction on,
a constant force of 8.0e-15 N per methane molecule along x, kerogen atoms
tethered by `spring/self` and held at 423 K by `langevin`, and the two
fixed barrier walls. The gas dump keeps its 11 columns
(`id type x y z vx vy vz c_dstress[1] c_dstress[2] c_dstress[3]`).
