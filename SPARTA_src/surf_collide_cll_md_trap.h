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

   cll/md/trap: cll/md with an explicit trapping-desorption branch. With
   probability p_trap(eps,theta), read from a second MD table, the
   particle is re-emitted thermally at Tsurf, otherwise it follows the
   cll/md path. The alpha table is calibrated on returned particles only,
   so it applies to the untrapped fraction alone.
   Usage: surf_collide <id> cll/md/trap Tsurf acc_rot acc_vib
          <alpha_table> <ptrap_table>
------------------------------------------------------------------------- */

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(cll/md/trap,SurfCollideCLLMDTrap)

#else

#ifndef SPARTA_SURF_COLLIDE_CLL_MD_TRAP_H
#define SPARTA_SURF_COLLIDE_CLL_MD_TRAP_H

#include "pointers.h"
#include "surf_collide.h"

namespace SPARTA_NS {

class SurfCollideCLLMDTrap : public SurfCollide {
 public:
  SurfCollideCLLMDTrap(class SPARTA *, int, char **);
  ~SurfCollideCLLMDTrap();
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

  // MD p_trap table, on its own grid (usually the same as the alpha grid)
  int pt_n_eps,pt_n_theta;
  double *pt_eps_grid;
  double *pt_theta_grid;
  double *ptrap_tbl;               // [pt_n_eps*pt_n_theta]

  double vstream[3];
  class RanKnuth *random;

  void cll_md_trap(Particle::OnePart *, double *);
  void cll_scatter(Particle::OnePart *, double *);
  void diffuse_emit(Particle::OnePart *, double *, double *, double *);
  void read_alpha_table(const char *fname);
  void read_ptrap_table(const char *fname);
  void interpolate_alpha(double eps_eV, double theta_deg,
                         double &a_t, double &a_n);
  double interpolate_ptrap(double eps_eV, double theta_deg);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.

E: Cannot open surf_collide cll/md/trap alpha table file

The named MD alpha table file could not be opened for reading.

E: Cannot open surf_collide cll/md/trap p_trap table file

The named MD p_trap table file could not be opened for reading.

*/
