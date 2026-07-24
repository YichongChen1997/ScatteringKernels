# SPARTA surface collision styles

Three `surf_collide` styles that take the reflection law from molecular dynamics
data instead of a fixed accommodation coefficient. They were written against the
SPARTA release of 24 Sep 2025.

## Build

Copy the six files into the SPARTA source directory and rebuild:

    cp surf_collide_cll_md.cpp surf_collide_cll_md.h \
       surf_collide_cll_md_trap.cpp surf_collide_cll_md_trap.h \
       surf_collide_dvsm_md.cpp surf_collide_dvsm_md.h  /path/to/sparta/src/
    cd /path/to/sparta/src
    make mpi

Each style registers itself through the usual `SurfCollideStyle` macro, so no
other source file has to be touched. `Tsurf` takes the same forms as in the
stock styles, including a per-surface temperature variable.

## cll/md

    surf_collide <id> cll/md Tsurf acc_rot acc_vib <alpha_table>

CLL reflection in which the normal and tangential accommodation coefficients are
looked up per collision from a table `alpha_t(eps,theta)`, `alpha_n(eps,theta)`,
indexed by the incident energy in eV and the angle to the local normal in
degrees. Rotational and vibrational accommodation stay constant and must lie in
[0,1]. The lookup is bilinear and clamped at the edges of the grid, so any
extrapolation past the MD points has to be built into the table itself. The
scattering is that of the stock `cll` style; its translate, rotate and partial
options are not supported.

## cll/md/trap

    surf_collide <id> cll/md/trap Tsurf acc_rot acc_vib <alpha_table> <ptrap_table>

The same as `cll/md`, with an explicit trapping and desorption branch in front of
it. With probability `p_trap(eps,theta)`, read from a second table on its own
grid, the particle is re-emitted thermally at `Tsurf` with full internal
accommodation; otherwise it follows the `cll/md` path. This split matters because
the alpha table is calibrated on returned particles only, so it describes the
untrapped fraction of the flux alone. A table of zeros reproduces `cll/md`.

## dvsm/md

    surf_collide <id> dvsm/md Tsurf <event_library>

Non-separable reflection with no kernel: outgoing velocities are resampled
directly from a library of MD beam events, in the spirit of the DVSM approach of
Mehta et al. (2017). The incident energy and angle select an event bin, the drawn
velocity is rotated into the local surface frame and rescaled by
`sqrt(eps/eps_bin)`, and the out-of-plane tangential sign is randomised. With the
per-bin probability `p_trap` the particle desorbs thermally at `Tsurf` instead.
The library carries no internal state, so rotational and vibrational energies
relax fully to `Tsurf`.

File formats for the two tables and the event library are given in
`../SPARTA_input/README.md`.
