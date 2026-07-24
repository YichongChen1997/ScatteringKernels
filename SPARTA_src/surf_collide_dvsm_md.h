/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Yichong Chen
   (Institute of Mechanics, Chinese Academy of Sciences)

   dvsm/md: non-separable, kernel-free surface collision that resamples
   outgoing velocities from a library of MD beam events (DVSM-style, cf.
   Mehta et al. 2017). The incident (eps,theta) selects an event bin, the
   drawn velocity is rotated into the local surface frame and rescaled by
   sqrt(eps/eps_bin), and with probability p_trap(bin) the particle
   desorbs thermally at Tsurf instead.
   Usage: surf_collide <id> dvsm/md Tsurf <event_library>
------------------------------------------------------------------------- */

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(dvsm/md,SurfCollideDVSMMD)

#else

#ifndef SPARTA_SURF_COLLIDE_DVSM_MD_H
#define SPARTA_SURF_COLLIDE_DVSM_MD_H

#include "pointers.h"
#include "surf_collide.h"

namespace SPARTA_NS {

class SurfCollideDVSMMD : public SurfCollide {
 public:
  SurfCollideDVSMMD(class SPARTA *, int, char **);
  ~SurfCollideDVSMMD();
  void init();
  Particle::OnePart *collide(Particle::OnePart *&, double &,
                             int, double *, int, int &);

 private:
  // MD event library on (eps,theta) grid
  int n_eps,n_theta;               // grid sizes
  double *eps_grid;                // [n_eps]  incident energy bins (eV)
  double *theta_grid;              // [n_theta] incidence angle bins (deg)
  int *nevents;                    // [nbin] outgoing events per bin
  double *ptrap;                   // [nbin] trapping probability per bin
  int *offset;                     // [nbin] start index into events (triplets)
  double *events;                  // [3*ntotal] vx,vy,vz per event (m/s),
                                   // beam frame: +z outward normal,
                                   // +x incident tangential

  double vstream[3];
  class RanKnuth *random;

  void dvsm(Particle::OnePart *, double *);
  void diffuse_emit(Particle::OnePart *, double *,
                    double *tangent1, double *tangent2);
  void read_event_library(const char *fname);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.

E: Cannot open surf_collide dvsm/md event library file

The named MD event library file could not be opened for reading.

*/
