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

   cll/md: CLL reflection in which acc_n and acc_t are looked up per
   collision from an MD table alpha_t(eps,theta), alpha_n(eps,theta),
   using the incident energy and the angle to the local normal.
   Usage: surf_collide <id> cll/md Tsurf acc_rot acc_vib <alpha_table>
------------------------------------------------------------------------- */

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(cll/md,SurfCollideCLLMD)

#else

#ifndef SPARTA_SURF_COLLIDE_CLL_MD_H
#define SPARTA_SURF_COLLIDE_CLL_MD_H

#include "pointers.h"
#include "surf_collide.h"

namespace SPARTA_NS {

class SurfCollideCLLMD : public SurfCollide {
 public:
  SurfCollideCLLMD(class SPARTA *, int, char **);
  ~SurfCollideCLLMD();
  void init();
  Particle::OnePart *collide(Particle::OnePart *&, double &,
                             int, double *, int, int &);

 private:
  double acc_rot,acc_vib;          // rot/vib accommodation, constant

  // MD alpha table: alpha_t(eps,theta), alpha_n(eps,theta)
  int n_eps,n_theta;               // grid sizes
  double *eps_grid;                // [n_eps]  incident energy bins (eV)
  double *theta_grid;              // [n_theta] incidence angle bins (deg)
  double *alpha_t_tbl;             // [n_eps*n_theta] flattened, ie*n_theta+ith
  double *alpha_n_tbl;             // [n_eps*n_theta]

  double vstream[3];
  class RanKnuth *random;

  void cll_md(Particle::OnePart *, double *);
  void read_alpha_table(const char *fname);
  void interpolate_alpha(double eps_eV, double theta_deg,
                         double &a_t, double &a_n);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.

E: Cannot open surf_collide cll/md alpha table file

The named MD alpha table file could not be opened for reading.

*/
