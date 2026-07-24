# Ar-Pt molecular beam runs

A monoenergetic argon beam is fired at a platinum slab and the outgoing
velocities are recorded. The runs are the calibration data behind the three
SPARTA wall models in `../../SPARTA_src`.

`in.beam` inserts argon atoms at a fixed speed and angle, lets them scatter, and
deletes them once they are back in the bulk. Gas-gas interaction is switched off,
so the atoms do not see each other and each one is an independent scattering
event. Interactions are Lennard-Jones with a 15 Angstrom cutoff. The wall is held
at `T_WALL` with a Nose-Hoover thermostat and is equilibrated for 10 ps before
the first insertion.

## Grid

Five incident energies, 0.4, 1, 2, 5 and 10 eV, at five polar angles from the
surface normal, 15, 30, 45, 60 and 75 degrees. That is 25 conditions, each with
2500 inserted atoms. The azimuth is fixed at PHI = 0, so the incident tangential
direction is +x.

`in.beam` takes the incident speed rather than the energy, through
`V_INC = sqrt(2 eps/m)`:

| eps [eV] | V_INC [m/s] |
|---------:|------------:|
|      0.4 |        1390 |
|      1   |        2198 |
|      2   |        3108 |
|      5   |        4915 |
|     10   |        6950 |

## Running

`in.beam` reads `../../ArPt_slab.data`, so it expects one directory per condition
under `runs/`. The slab file is in this directory: 169744 atoms, type 1 for the
gas, 2 for the thermostatted wall and 3 for the fixed layer. A different wall can
be built with the tools in `../../../initialisation/TypeOfWalls`.

Name each directory `runs/eps<E>_th<theta>`, since the post-processing reads the
energy and angle back out of the directory name:

    mkdir -p runs/eps10_th45 && cd runs/eps10_th45
    mpirun -np 4 lmp_mpi -in ../../in.beam \
        -var V_INC 6950 -var THETA 45 -var PHI 0 \
        -var T_WALL 300 -var T_GAS 300 \
        -var n_insert 2500 -var n_loops 1

The product of `n_insert` and `n_loops` has to be 2500, because the trapping
probability is formed as `1 - N_return/2500` downstream. `TAILSTEPS`, 8000 steps
by default, is the run after the last insertion, in which the atoms still in
flight finish scattering and reach the dump. It has to be raised at low incident
energy, where the atoms are slow. A finished run prints `Simulation completed`.

The default `velocityMode 2` gives the monoenergetic beam used here.
`velocityMode 1` adds a thermal spread at `T_GAS`. The default `wallMode 1`
thermostats the whole wall; `wallMode 2` leaves a free surface layer of thickness
`ablationThick` on a thermostatted base.

## From runs to tables

    in.beam, 25 conditions
      -> runs/eps*_th*/dump_gas.lammpstrj
    extract_alpha.py --all
      -> calibration_curve_ArPt.csv, runs/eps*_th*/events_out.csv
    build_table_arpt.py --install
      -> ../../SPARTA_input/tables/alpha_arpt_table{,_hold}.txt
    build_ptrap_table.py --install
      -> ../../SPARTA_input/tables/ptrap_arpt_table{,_hold}.txt
    build_event_library.py --install
      -> ../../SPARTA_input/events/dvsm_events.txt

`extract_alpha.py` takes an atom as returned if its last appearance in the dump
has z above the 15 Angstrom cutoff and vz > 0, where the wall exerts no force and
the velocity is the free-flight one. The tangential and normal coefficients are
averaged over the returned atoms only, with a bootstrap standard error, and the
count of returned atoms is carried in the `events` column. The two `build_*`
scripts and the event library all read that column, so the coefficients, the
trapping probability and the event bins stay consistent.

Without `--install` the scripts write the tables here instead of into the SPARTA
directories.

## Caveats in the tables

The measured 5 by 5 grid enters the tables unchanged, so no separable form
`alpha(eps,theta) = g(eps) h(theta)` is assumed. It is then padded: 0.05 eV holds
the 0.4 eV row, 0 degrees holds the 15 degree column, and 90 degrees is either a
linear extrapolation of the 60 and 75 degree values or a hold of the 75 degree
value. The two variants bracket that extrapolation. Energies above 10 eV are held
at the 10 eV values by the lookup itself.

The CL kernel needs coefficients in [0,1], so the tables are clipped. The normal
coefficient comes out negative at some grazing conditions and is clipped to zero
there, which biases normal accommodation upward. Clipped points are listed when
the table is built.

At 75 degrees and 0.4, 1, 2 and 5 eV no atom returns at all. The normal energy is
too low for anything to escape the well, and there is nothing to average over.
Those four points are set to the complete accommodation limit rather than
measured, on the grounds that trapped atoms desorb thermalised. They carry
p_trap = 1 in the trapping table, so `cll/md/trap` and `dvsm/md` never reach the
coefficients there.
