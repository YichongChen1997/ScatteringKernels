# SPARTA cases

Two 2D argon cases at 7 km/s: flow over a cylinder (`in.cyl`, D = 0.1 m) and over
a flat plate at 20 degrees incidence (`in.plate`, chord 0.1 m, thickness 6 mm).
Both read the wall models built from `../SPARTA_src`.

## Running a case

Write the surface geometry, write the run directories, then execute one:

    python3 make_geom.py
    python3 gen_runs.py
    runs/cyl_kn1_cllmd_extrap/cmd.sh

Paths inside the decks are relative to the run directory, two levels below this
one, so run a deck from `runs/<case>/` rather than from here. `cmd.sh` holds the
full command line if you would rather run it by hand.

`make_geom.py` writes `geom/circle_D0.1.surf` and `geom/plate_c0.1_a20.surf`.
Points run clockwise. With the opposite ordering SPARTA treats the inside of the
body as the flow region and emits no particles.

Variables that can be set on the command line: `SEED`, `WALLMODEL`, `SIGMA`,
`ALPHAT`, `ALPHAN`, `ALPHATABLE`, `PTRAPTABLE`, `NRHO`, `UINF`, `TINF`, `TWALL`,
`NPC`, `NX`, `NY`, `DT`, `NUNSTEADY`, `NSTEADY`. `NUNSTEADY` and `NSTEADY` must
be multiples of 100, because the stats output reduces a fix with Nfreq = 100.

`gen_runs.py` writes one directory per case under `runs/`, covering 2 geometries
x 3 Knudsen numbers (0.1, 1, 10) x 7 wall-model settings, each with a `cmd.sh`,
and a `manifest.tsv`. With `--seeds 4` it repeats every case under `runs_ver/`
with different random seeds, which is what the confidence intervals below need. Kn is converted to a free-stream number density through the VHS
mean free path with the argon parameters of `ar.vss`, taking D as the cylinder
diameter or the plate chord. Set `SPARTA_BIN` if the executable is not called
`spa_mpi`.

Each run writes `surf_force.dump`, holding per-facet `fx`, `fy` in N per metre of
span and the net energy flux `etot` in W/m^2.

## Wall models

`-var WALLMODEL` selects one of five boundary conditions, from a fixed reflection
law to direct use of the MD data:

| WALLMODEL   | boundary condition                                     | extra variables        |
|-------------|--------------------------------------------------------|------------------------|
| `diffuse`   | diffuse at `TWALL` with accommodation `SIGMA`          | `SIGMA`                |
| `cll`       | CLL with constant coefficients                          | `ALPHAN`, `ALPHAT`     |
| `cllmd`     | CLL with `alpha_n`, `alpha_t` resolved in energy and angle | `ALPHATABLE`        |
| `cllmdtrap` | as `cllmd`, plus an explicit trapping branch            | `ALPHATABLE`, `PTRAPTABLE` |
| `dvsm`      | direct resampling of the MD event library               | `ALPHATABLE`           |

`dvsm` passes the event library through `ALPHATABLE`, so use
`-var ALPHATABLE events/dvsm_events.txt`.

The supplied Ar-Pt data are in `tables/` and `events/`. Each table comes in two
variants, `_table.txt` and `_table_hold.txt`, which differ only at grazing
incidence: one extrapolates the 60 and 75 degree values to 90 degrees, the other
holds the 75 degree value. Running both brackets that extrapolation.

## Post-processing

Run these from this directory.

`tools/analyze.py` reads `runs/*/surf_force.dump` and reports C_D and the
stagnation heat flux for the cylinder, and C_L, C_D and C_m about the leading
edge for the plate. `--all` scans `runs/` and writes `runs/summary.tsv`.

`tools/analyze_ci.py` combines four seeds per case into a 95 per cent confidence
interval and pairs the deviation against the diffuse baseline by seed. It also
reports a refinement check in timestep and particles per cell, and a degeneracy
check that `cll/md` with a constant table reproduces the native `cll` statistics.
The replicates come from `gen_runs.py --seeds 4`. The refined and degeneracy runs
are not generated: copy a `cmd.sh` into `runs_ver/conv_<geom>_kn1_cllmd_dt2`,
`..._npc2` and `runs_ver/degen_<geom>_kn1_cllmdconst`, and add `-var DT 5e-8`,
`-var NPC 40` or a constant accommodation table respectively. Writes
`runs_ver/ci_summary.tsv`.

`tools/analyze_dvsm.py` separates the effect of the energy-resolved coefficients
from the shape of the kernel, by comparing the separable CLL fit against the
resampled MD events. Writes `runs_ver/dvsm_decomposition.tsv`.

`tools/fm_analytic.py` gives free-molecular coefficients from the panel formulae
of Schaaf and Chambre (1961), as an analytic check on the Kn = 10 runs.

## File formats

Comment lines start with `#`. Values are read as a plain token stream, so line
breaks carry no meaning. Both grids must be ascending and hold at least two
points.

Accommodation table, read by `cll/md` and `cll/md/trap`:

    n_eps n_theta
    eps_grid       n_eps values, eV
    theta_grid     n_theta values, deg
    alpha_t        n_eps*n_theta values, eps outer, theta inner
    alpha_n        n_eps*n_theta values

Trapping table, read by `cll/md/trap`, on its own grid, which need not match the
alpha grid:

    n_eps n_theta
    eps_grid       n_eps values, eV
    theta_grid     n_theta values, deg
    ptrap          n_eps*n_theta values, all in [0,1]

Event library, read by `dvsm/md`:

    n_eps n_theta
    eps_grid       n_eps values, eV
    theta_grid     n_theta values, deg
    per bin, eps outer and theta inner:
      eps theta N_return p_trap
      N_return lines of  vx vy vz

Event velocities are in m/s in the beam frame, with +z along the outward normal
and +x along the incident tangential direction. Bin headers must appear in grid
order and match the declared grids, and a bin with no events must have
p_trap = 1. The supplied library has 25 bins and 50251 events.

For another gas-surface pair, supply tables on the same format. Lookup in
`cll/md` is clamped at the grid edges and the coefficients are clipped to [0,1],
so put whatever behaviour you want outside the measured range into the end
points. Event velocities carry the mass of the species they were generated with,
so a library is tied to its gas.
