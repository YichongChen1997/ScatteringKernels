# The post-processing tools in `tools/`

`tools/` holds 21 small C++ programs, `pp_meas_*.cpp`, that read the dump
files of a LAMMPS run of gas between walls and write accommodation
coefficients, velocity and angle distributions, density profiles and a few
other things. They were used for the published papers and the thesis, so
their source is frozen: it is not changed, bugs included. This page says
what each program does, what it needs, what it writes, and what is known to
be wrong with it.

Line numbers below (for example `pp_meas_ACs.cpp:184`) refer to the files in
`tools/`.

## Building and running

```bash
make                       # compiles every tools/*.cpp into build/bin/
```

Each program is compiled on its own with `g++ -o name name.cpp` and no
flags, which is how the published results were made. Use `make CXX=clang++`
or similar to pick another compiler.

To run a tool, go to a directory that holds the dump file(s) and a
`Specification.dat`, and start it from there:

```bash
cd my_run                  # holds dump_meas_gas.lammpstrj (and dump_meas_wall.lammpstrj)
cp /path/to/ScatteringKernels/tools/Specification.example.dat Specification.dat
# edit Specification.dat, then put the number of frames on its first line:
make -f /path/to/ScatteringKernels/scripts/rundir.mk update
/path/to/ScatteringKernels/build/bin/pp_meas_ACs
```

The programs take no arguments and read nothing from the keyboard. They
write their output files into the current directory, overwriting files of
the same name, and print their progress (one or two lines per frame) and a
short summary to the screen.

`make test` runs every tool on the small stored reference sets in
`tests/reference/` and compares the output with the stored answers (see
`tests/reference/README.md`).

## What the tools expect

### Units

The dump must come from a run with LAMMPS `units real`: lengths in Å, time
in fs, velocities in Å/fs, energies in kcal/mol, forces in kcal/mol/Å. The
factors are written into each program (for example 1 Å/fs = 10^5 m/s,
`pp_meas_ACs.cpp:22-31`). A dump in any other unit system is read without
complaint and gives wrong numbers.

### The dump files

- The gas dump must be called `dump_meas_gas.lammpstrj` and the wall dump
  (only `pp_meas_wallTemp` reads one) `dump_meas_wall.lammpstrj`.
- Every frame has the usual 9 header lines. The tools take the time step
  from line 2, the number of atoms from line 4 and the box from lines 6-8;
  they do not look at the `ITEM: ATOMS` line.
- The gas dump must have exactly these 16 columns, in this order:
  `id type x y z vx vy vz c_dstress[1] c_dstress[2] c_dstress[3] c_KE c_PE c_FORCE[1] c_FORCE[2] c_FORCE[3]`,
  as written by `LAMMPS_input/Explicit_layered/in.meas` (line 100). The
  three stress columns are the diagonal of the per-atom stress (compute
  `stress/atom`, pressure times volume), `c_KE` and `c_PE` the per-atom
  kinetic and potential energy, `c_FORCE` the per-atom force. See KI-8 for
  what happens with other columns.
- The wall dump must have 8 columns: `id type x y z vx vy vz`.
- Atoms must come in the same order in every frame and their number must
  not change (`dump_modify sort id` in the deck). The collision tools keep
  track of an atom by its line number within the frame, not by its id.

### Gas atoms are type 1

The tools that follow molecules count only atoms of type 1 as gas
(for example `pp_meas_ACs.cpp:186`). `pp_meas_Bins` and
`pp_meas_GasOverTime` use every atom in the gas dump, whatever its type, so
the gas dump should contain gas only.

### The virtual plane at z = rCut

The bottom wall's surface is taken to be at z = 0 (as in
`initialisation/TypeOfWalls/Explicit_layered`). A virtual plane at
z = rCut (line 3 of `Specification.dat`) separates the gas from the
near-wall layer. A gas-wall "collision" is one visit of a molecule below
this plane (`pp_meas_ACs.cpp:186-277`, and the same code in the other
collision tools):

1. A molecule counts only once it has been seen above the plane
   (z >= rCut), in the first frame or later.
2. When it is first seen below the plane (z < rCut), the collision has
   started. Its **incident** position and velocity are those of the last
   frame in which it was still above the plane.
3. When it is next seen above the plane again, the collision has ended.
   Its **reflected** position and velocity are those of this first frame
   back above the plane.

So a collision is only counted once it has ended; molecules still below the
plane in the last frame are dropped (this is why "Collisions Started" is
larger than "Collisions Ended" on the screen). Only the bottom wall is
analysed: "No. of Collisions at top" is always 0. The plane also counts
gas-gas collisions below z = rCut as wall collisions if gas atoms interact
with each other (see KI-15).

In the output, "normal" means |vz|, "tangential x" means vx (with its
sign), and "tangential" alone means sqrt(vx² + vy²). Many results are
grouped by the incident speed in units of

    vM = sqrt(2 kB Tw / mG),

the most probable speed of the gas at the wall temperature Tw (about
353 m/s for Ar at 300 K).

### Specification.dat

Five lines, read by position. On each line the tools read the numbers they
need and skip the rest of the line, so anything after the numbers (such as
a `# comment`) is ignored. There are no names; the order is what counts.
`tools/Specification.example.dat` is a filled-in example.

| Line | Values | Meaning | Used by |
|---|---|---|---|
| 1 | `nTimeSteps` | number of frames to read from the gas dump; must equal the number of frames in the file (KI-16). `make -f .../scripts/rundir.mk update` in the run directory writes it for you | all |
| 2 | `deltaT tSkip skipTimeStep` | MD time step in fs; MD steps between two frames of the gas dump; frames left out at the start | `deltaT` and `tSkip` only to label time in `pp_meas_GasOverTime`, `pp_meas_ResidenceTime`, `pp_meas_wallTemp`; `skipTimeStep` only in `pp_meas_Bins` and `pp_meas_GasOverTime` |
| 3 | `rCut H` | height of the virtual plane above the bottom surface, Å; channel height (position of the top surface), Å | `rCut` by all collision tools, `pp_meas_Contours`, `pp_meas_NearWall`; `H` by `pp_meas_GasOverTime` and `pp_meas_NearWall` |
| 4 | `Tg Tw` | gas and wall temperature, K | only `Tw` is used (for vM); `Tg` is read but never used |
| 5 | `mG mW` | mass of one gas atom and one wall atom, g/mol | `mG` by most tools; `mW` only by `pp_meas_wallTemp` |

If `Specification.dat` is missing the tools stop with "Error opening
Specification.dat file". A missing dump file is **not** reported: the tools
run to the end and write empty or meaningless files.

### Histograms

Most distributions are written with three columns: bin centre, number of
events in the bin, and that number divided by the area under the histogram
(trapezoid rule), so that the third column integrates to about 1. Velocity
histograms use 10 m/s bins up to 4 vM (142 bins for Ar at 300 K); a signed
quantity such as vx gets the negative half first, then the positive half.
Where a file holds several distributions (one per window), they are written
one after the other with no blank line between them.

## The tools at a glance

| Tool | What it gives | Writes | Known issues |
|---|---|---|---|
| `pp_meas_ACs` | accommodation coefficients (TMAC, NMAC, tangential and normal energy, total energy), two formulas, all collisions and by incident-speed window | `AlphaTx.txt` `AlphaN.txt` `AlphaEx.txt` `AlphaEn.txt` `AlphaE.txt` | KI-8, KI-15, KI-16 |
| `pp_meas_AngularF` | incident and reflected polar angle, all collisions | `thetaPrimeDis_overall.txt` `thetaDis_overall.txt` | KI-8, KI-16 |
| `pp_meas_AngularP` | reflected polar angle for five incident-angle windows | `thetaDis_beam.txt` | KI-8, KI-16 |
| `pp_meas_Bins` | profiles across the channel: density, velocity, temperature, pressure | `bins_meas_coarse.txt` `bins_meas_fine.txt` | KI-8, KI-16 |
| `pp_meas_Collisions` | how many times a molecule bounces while below the plane, and TMAC / NEAC by that number | `collisionDis.txt` `collisionAlphaTx.txt` `collisionAlphaEn.txt` | KI-8, KI-16 |
| `pp_meas_Contours` | 3-D maps of the gas potential energy and normal force near the bottom wall | `PES.txt` `Fz.txt` (300–400 MB each) | KI-8, KI-16 |
| `pp_meas_Correlation` | incident and reflected velocity of every collision | `bottomVelocities.txt` | KI-8, KI-16 |
| `pp_meas_Displacement` | distance travelled along the wall during a collision | `displacementDis.txt` | KI-5, KI-8, KI-16 |
| `pp_meas_FAngular` | same as `pp_meas_AngularF`, other file names | `Incident_angular_full.txt` `Reflected_angular_full.txt` | KI-8, KI-16 |
| `pp_meas_GasOverTime` | gas velocity, flow rate, temperature and pressure against time | `gasOverTime.txt` | KI-4, KI-8, KI-16 |
| `pp_meas_NearWall` | positions of gas atoms near either wall at four moments | `nearWallMolecules.xyz` | KI-8, KI-16 |
| `pp_meas_PAngular` | same as `pp_meas_AngularP`, other file name | `reflected_angle_partial.txt` | KI-8, KI-16 |
| `pp_meas_ResidenceTime` | time spent below the plane per collision | `timeDis.txt` | KI-8, KI-16 |
| `pp_meas_Traj` | trajectories of four chosen atoms | `molecules.trj` `distanceToSurf.txt` | KI-7, KI-8, KI-16 |
| `pp_meas_VDFF` | incident and reflected velocity distributions, all collisions | `xi_Nprime_overall.txt` `xi_N_overall.txt` `xi_Tprime_overall.txt` `xi_T_overall.txt` `xi_MagPrime_overall.txt` `xi_Mag_overall.txt` `xi_TxPrime_overall.txt` `xi_Tx_overall.txt` | KI-1, KI-8, KI-16 |
| `pp_meas_VDFPi` | incident velocity distributions by incident-speed window | `xi_Nprime_beam.txt` `xi_TxPrime_beam.txt` | KI-2, KI-8, KI-16 |
| `pp_meas_VDFPr` | reflected velocity distributions by incident-speed window | `reflected_N.txt` `reflected_Tx.txt` | KI-8, KI-16 |
| `pp_meas_penACs` | TMAC and NEAC by how deep the molecule got | `depthAlphaTx.txt` `depthAlphaEn.txt` | KI-3, KI-8, KI-16 |
| `pp_meas_penBeam` | reflected velocity distributions by depth and incident-speed window | `depthBeamDis_Tx.txt` `depthBeamDis_N.txt` | KI-8, KI-16 |
| `pp_meas_penOverall` | depth distribution, and reflected velocity distributions by depth | `depthDis.txt` `depthOverallDis_Tx.txt` `depthOverallDis_N.txt` | KI-8, KI-16 |
| `pp_meas_wallTemp` | temperature of the wall layers against time, from the wall dump | `MeasTemp_BottomWall.txt` `MeasTemp_TopWall.txt` | KI-6, KI-8, KI-16 |

All read `Specification.dat` and `dump_meas_gas.lammpstrj`, except
`pp_meas_wallTemp`, which reads `Specification.dat` and
`dump_meas_wall.lammpstrj`. The file names above are the ones in the code;
`tools/README.txt` (the author's original notes, kept as they were) has a
few that are out of date, for example `xi_N_beam.txt` for what
`pp_meas_VDFPr` writes as `reflected_N.txt`.

## The tools one by one

### pp_meas_ACs

Accommodation coefficients of the bottom wall from the incident and
reflected velocities of every collision: TMAC from vx (`AlphaTx`), NMAC
from |vz| (`AlphaN`), tangential-x energy from vx² (`AlphaEx`), normal
energy from vz² (`AlphaEn`) and total energy from |v|² (`AlphaE`). Each is
worked out in two ways: the "general formula" (ratio of summed incident
and reflected values, thesis eq. 3.6) and a least-squares fit of reflected
against incident values (thesis eq. 3.7).

Each file has 29 rows, one per window of incident speed centred on
0.1, 0.2, ..., 2.9 vM, ±0.1 vM (for `AlphaE` the window is on |v|² in
units of kB Tw / mG). Columns:

1. window centre;
2. general formula, collisions in this window;
3. least squares, collisions in this window;
4. general formula, all collisions;
5. least squares, all collisions.

Columns 4 and 5 are the same on every row. The general formula for TMAC
uses only collisions with incident vx >= 0, and for NMAC only those with
|vz| < 6 vM (`pp_meas_ACs.cpp:367-379`). Empty windows give `nan`. The
screen shows how many collisions started and ended.

### pp_meas_AngularF and pp_meas_FAngular

Distribution of the polar angle from the surface normal,
atan(sqrt(vx² + vy²) / |vz|), in 1° bins (180 rows; only 0-90° can be
filled). `pp_meas_AngularF` writes the incident angle to
`thetaPrimeDis_overall.txt` and the reflected angle to
`thetaDis_overall.txt`; `pp_meas_FAngular` is the same program writing
`Incident_angular_full.txt` and `Reflected_angular_full.txt`
(`pp_meas_AngularF.cpp:381,428`, `pp_meas_FAngular.cpp:381,428`; the two
files differ only in these names).

### pp_meas_AngularP and pp_meas_PAngular

Distribution of the reflected polar angle for collisions whose incident
angle is within 2° of 15°, 30°, 45°, 60° and 75°: five histograms of 180
rows, one after the other, in `thetaDis_beam.txt` (`pp_meas_AngularP`) or
`reflected_angle_partial.txt` (`pp_meas_PAngular`, otherwise identical).
The screen shows how many collisions fall in each window.

### pp_meas_Bins

Profiles along z over the whole box height, averaged over all frames from
`skipTimeStep` on, in 0.5 Å bins (`bins_meas_coarse.txt`) and 0.25 Å bins
(`bins_meas_fine.txt`); the bin width is adjusted slightly so that the box
holds an even number of bins. Columns: z of the bin centre (Å), mass
density (kg/m³), number density (1/m³), mean vx (m/s), temperature (K,
from the velocity minus the bin's mean vx), pressure (MPa, from the per-atom
stress). Every atom in the dump is used, whatever its type.

### pp_meas_Collisions

While a molecule is below the plane, the tool counts how often vz changes
sign (a "bounce"). `collisionDis.txt`: histogram of the number of bounces
per collision, 0 to 19 (bin, count, normalised). `collisionAlphaTx.txt`:
for 0 to 19 bounces, TMAC by least squares (column 2) and by the general
formula (column 3). `collisionAlphaEn.txt`: normal energy accommodation by
the general formula (column 2) and least squares (column 3). Note that the
two files have the two methods in opposite order
(`pp_meas_Collisions.cpp:450-454, 485-498`). Rows with no collisions give
`nan`.

### pp_meas_Contours

Averages the potential energy (converted to K per atom) and the z force (N)
of gas atoms below the plane on an 800 × 800 × (20 × rCut) grid over x, y
and z, and writes `PES.txt` and `Fz.txt`: for each z layer, 800 lines of
800 numbers, then a blank line. Empty cells are 0. The grid grows with
rCut: at rCut = 12 Å (R3) it needs about 4 GB of memory and writes about
300 MB per file, at rCut = 15 Å (R2) about 5.5 GB and 384 MB, whatever the
number of frames.
The x and y cells assume the box starts at 0. `make test` leaves it out
unless asked (see the Makefile).

### pp_meas_Correlation

One line per collision in `bottomVelocities.txt`: incident vx, vy, vz, then
reflected vx, vy, vz, all in m/s. This is the raw data behind the scatter
plots of reflected against incident velocity (thesis Fig. 3.4).

### pp_meas_Displacement

Distance in the x-y plane between where a molecule was when the collision
started (last frame above the plane) and where it was when it ended (first
frame back above), as a histogram in 0.5 Å bins up to 200 Å
(`displacementDis.txt`). See KI-5 for the periodic boundaries.

### pp_meas_GasOverTime

One line per frame in `gasOverTime.txt`: time since the first frame (s),
mean vx of the gas (m/s), mass flow rate along x (kg/s), the sum of
m |v - mean vx|² over all atoms in kg Å²/fs² (twice the thermal kinetic
energy; multiply by 10^10 for J), temperature (K), and pressure (MPa). The
pressure is a running average over the frames since `skipTimeStep`, using
the volume Lx × Ly × H. Every atom in the dump is used, whatever its type.
The pressure column is wrong because of KI-4.

### pp_meas_NearWall

Positions of the gas atoms that are within rCut of either wall
(z < rCut or z > H - rCut) in four frames: numbers 2N/3, N/2, 2N/5 and N/3
(integer division, counting from 0), where N is `nTimeSteps`.
`nearWallMolecules.xyz`: the total count, the count in each of the four
frames, then one line per atom: frame index (0-3), x, y, z (Å).

### pp_meas_ResidenceTime

How long each collision lasted, from the last frame above the plane to the
first frame back above, as a histogram in steps of one frame, up to 200
frames (`timeDis.txt`). Columns: bin centre in s, count, normalised, and
the mean duration in s (the same on every row).

### pp_meas_Traj

Follows the atoms with ids 53, 235, 1000 and 5000, fixed in the code
(`pp_meas_Traj.cpp:57-60`), and writes their positions as a dump that OVITO
can read (`molecules.trj`) and their z against frame number
(`distanceToSurf.txt`). If any of these ids is missing from the dump it
crashes (KI-7).

### pp_meas_VDFF

Velocity distributions over all collisions, incident (files with "prime")
and reflected: normal |vz| (`xi_Nprime_overall`, `xi_N_overall`),
tangential sqrt(vx² + vy²) (`xi_Tprime_overall`, `xi_T_overall`), speed
(`xi_MagPrime_overall`, `xi_Mag_overall`, wrong, see KI-1) and signed vx
(`xi_TxPrime_overall`, `xi_Tx_overall`, twice as many rows). The screen
shows the mean incident and reflected vx.

### pp_meas_VDFPi

Incident velocity distributions for 20 windows of incident speed centred on
0.1, 0.2, ..., 2.0 vM (±0.1 vM): |vz| selected by incident |vz|
(`xi_Nprime_beam.txt`, 20 histograms) and vx selected by incident vx
(`xi_TxPrime_beam.txt`, 20 signed histograms). The second file is wrong,
see KI-2.

### pp_meas_VDFPr

The same windows as `pp_meas_VDFPi`, but the histograms are of the
**reflected** velocity: reflected |vz| for each window of incident |vz|
(`reflected_N.txt`), and reflected vx for each window of incident vx
(`reflected_Tx.txt`). The tangential part is written correctly here.

### pp_meas_penACs

For every collision the tool keeps the lowest z the molecule reached while
below the plane (its "penetration depth"; negative means below the surface
atoms' z = 0). It then works out TMAC and the normal energy accommodation
for depth windows about 2 Å wide covering -rCut to rCut (for rCut = 12 Å,
12 windows centred on -11, -9, ..., 11) and writes `depthAlphaTx.txt` and
`depthAlphaEn.txt` (columns: window centre, general formula, least
squares). The least-squares column is wrong, see KI-3. The
depth is updated one frame late (`pp_meas_penACs.cpp:231-234` compares with
the previous frame's z), so the last frame below the plane is not included.

### pp_meas_penBeam

Reflected velocity distributions for 10 depth windows (centres evenly
spaced from -rCut to rCut, each ±rCut/10 wide, so they do not quite touch)
times 10 windows of incident speed (0.2, 0.4, ..., 2.0 vM, ±0.2 vM, so
they overlap): reflected vx by incident vx (`depthBeamDis_Tx.txt`) and
reflected |vz| by incident |vz| (`depthBeamDis_N.txt`). The histograms are
written depth window by depth window, and within each depth window by
speed window.

### pp_meas_penOverall

`depthDis.txt`: histogram of the penetration depth (as in `pp_meas_penACs`)
in 1 Å bins from -rCut to rCut. `depthOverallDis_Tx.txt` and
`depthOverallDis_N.txt`: reflected vx and |vz| distributions for the same
10 depth windows as `pp_meas_penBeam`, all incident speeds together.

### pp_meas_wallTemp

Reads the wall dump and writes, for each frame after the first, the time
and the temperature of the surface layer (types 5 bottom / 6 top), of the
thermostatted bulk layers (types 2 / 3) and of both together, for the
bottom (`MeasTemp_BottomWall.txt`) and top (`MeasTemp_TopWall.txt`) wall.
Frozen atoms (type 4) are left out. The temperature is m v² / (3 kB)
averaged over the atoms, without removing any motion of the layer as a
whole. It reads only `nTimeSteps / 100` frames (KI-6).

## Known issues

These are recorded, not fixed: the published results were made with this
code. The numbers are the ones used for the project's list of known issues;
the numbers not used here belong to other parts of the repository.

**KI-1. `pp_meas_VDFF`: the speed |v| is wrong.** At
`pp_meas_VDFF.cpp:296-301`, `vMagT` is already vx² + vy², and the speed is
computed as sqrt(vz² + vMagT²), so the tangential part enters as
(vx² + vy²)² instead of vx² + vy². `xi_MagPrime_overall.txt` and
`xi_Mag_overall.txt` are wrong; the other six files are not affected.

**KI-2. `pp_meas_VDFPi`: the tangential file is mixed up.** For every
window two histograms (negative and positive vx) are appended to one list
(`pp_meas_VDFPi.cpp:461-464`), but the output reads entry `n` of that list
for window `n`, for both halves (`pp_meas_VDFPi.cpp:468,473`). Window `n`
therefore shows one half of window n/2, mirrored on both sides. Since the
windows only select positive incident vx, the negative halves are empty:
every even-numbered block of `xi_TxPrime_beam.txt` is all zeros, the odd
ones repeat the positive half of windows 0-9, and windows 10-19 never
appear. `xi_Nprime_beam.txt` is fine, and `pp_meas_VDFPr` has the fixed
version of this code.

**KI-3. `pp_meas_penACs`: the least-squares values are not per depth.**
The means are taken over the collisions in the depth window
(`pp_meas_penACs.cpp:355-380`), but the least-squares sums run over all
collisions (`pp_meas_penACs.cpp:391-395` and `409-413`). Column 3 of
`depthAlphaTx.txt` and `depthAlphaEn.txt` is therefore not the value for
that depth.

**KI-4. `pp_meas_GasOverTime`: `stressTemp` is never set to zero.** It is
declared at `pp_meas_GasOverTime.cpp:21` and added to at line 124 without
a starting value, so the pressure (column 6 of `gasOverTime.txt`, line 146)
starts from whatever was in memory. In the stored reference answers (Apple
clang) it happened to start at zero; built with GCC 15 the whole column is
`nan`. `make test` only warns about this column.

**KI-5. `pp_meas_Displacement`: periodic boundaries handled on one side
only.** `pp_meas_Displacement.cpp:295-303` corrects a jump in x or y larger
than half the box, but not one smaller than minus half the box. A molecule
that crosses the boundary in the negative direction gets a displacement of
nearly a box length.

**KI-6. `pp_meas_wallTemp`: the wall dump is assumed to be 100 times
sparser than the gas dump.** `pp_meas_wallTemp.cpp:60` divides
`nTimeSteps` by 100, and lines 145 and 151 label time with the gas dump's
`tSkip`. If the wall dump is written at any other interval, the tool reads
too few frames or reads past the end of the file (without an error), and
the time column is wrong. In `tests/reference/R3_channel_mini` the wall dump
has one frame, so the output is rows of zeros (full run) or empty (80
frames).

**KI-7. `pp_meas_Traj`: crashes when a tracked id is missing.** The ids 53,
235, 1000 and 5000 are fixed in the code (`pp_meas_Traj.cpp:57-60`). If one
of them is not a gas atom in the dump, its list stays empty and the output
loop reads past its end (`pp_meas_Traj.cpp:144`, and 161). The tool writes a
partial `molecules.trj` and then stops with a segmentation fault (exit
code 139), or, with a compiler that checks vector bounds (GCC 15 without
optimisation), with an abort (exit code 134).

**KI-8. Only the 16-column gas dump is read correctly, and nothing checks
it.** Every tool reads exactly 16 numbers per atom (for example
`pp_meas_ACs.cpp:184`; `pp_meas_wallTemp.cpp:91` reads 8 from the wall
dump) and never looks at the `ITEM: ATOMS` line. A dump with other columns,
such as the 11-column dumps of `examples/Fixed_adsorbate`, is misread
without any message: the reading falls out of step within the first frame,
after which every frame repeats the same time step and the results are
zeros or `nan`, and the tool still exits normally.

**KI-15. `LAMMPS_input/Explicit_layered/in.meas`: gas atoms interact with
each other.** Lines 51-52 of this deck switch the Ar-Ar interaction on (the
line that switches it off is there, commented out). With a mean free path
of about 4 Å, shorter than the 12 Å virtual plane, many of the
"collisions" the tools count below the plane are Ar-Ar collisions, and all
four accommodation coefficients come out near 0.95 in the reference case.
The thesis workflow switches gas-gas interaction off;
`tests/reference/R3_channel_mini/EDITS.md` shows the edit. The frozen deck
is left as it is.

**KI-16. Line 1 of `Specification.dat` must equal the number of frames.**
Each tool reads exactly `nTimeSteps` frames (for example
`pp_meas_ACs.cpp:50` and `121`) and never checks for the end of the file.
If the number is too small, the last frames are ignored. If it is too
large, the tool keeps "reading" past the end, reuses the last values it
read, and still exits normally. On the 80-frame reference dump with
`nTimeSteps = 100`, `pp_meas_ACs` prints time step 65800 twenty extra times
and reports 298 collisions started instead of 137. Run
`make -f /path/to/ScatteringKernels/scripts/rundir.mk update` in the run
directory to write the right number.
