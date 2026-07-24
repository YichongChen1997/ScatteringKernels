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

   cll/md: CLL reflection with acc_n and acc_t taken per collision from an
   MD table alpha(eps,theta). The scattering is that of the stock cll
   style; the translate/rotate/partial options of cll are not supported.
   Usage: surf_collide <id> cll/md Tsurf acc_rot acc_vib <alpha_table>
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "surf_collide_cll_md.h"
#include "surf.h"
#include "surf_react.h"
#include "particle.h"
#include "collide.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "math_extra.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,DISCRETE,SMOOTH};       // matches collide rotstyle/vibstyle enums

#define EV2J 1.602176634e-19      // 1 eV in J
#define MAXLINE 8192

/* ---------------------------------------------------------------------- */

SurfCollideCLLMD::SurfCollideCLLMD(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  // surf_collide <id> cll/md Tsurf acc_rot acc_vib <alpha_table>

  if (narg != 6) error->all(FLERR,"Illegal surf_collide cll/md command");

  parse_tsurf(arg[2]);

  acc_rot = atof(arg[3]);
  acc_vib = atof(arg[4]);

  if (acc_rot < 0.0 || acc_rot > 1.0 || acc_vib < 0.0 || acc_vib > 1.0)
    error->all(FLERR,"Surf_collide cll/md acc_rot/acc_vib must be >= 0 and <= 1");

  eps_grid = theta_grid = alpha_t_tbl = alpha_n_tbl = NULL;
  n_eps = n_theta = 0;
  read_alpha_table(arg[5]);

  vstream[0] = vstream[1] = vstream[2] = 0.0;

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideCLLMD::~SurfCollideCLLMD()
{
  delete random;
  delete [] eps_grid;
  delete [] theta_grid;
  delete [] alpha_t_tbl;
  delete [] alpha_n_tbl;
}

/* ---------------------------------------------------------------------- */

void SurfCollideCLLMD::init()
{
  SurfCollide::init();
  check_tsurf();
}

/* ----------------------------------------------------------------------
   read the MD alpha table on rank 0, broadcast to all ranks.
   File format (# starts a comment, tokens may be split over lines):
     n_eps n_theta
     eps_grid (n_eps floats, eV)
     theta_grid (n_theta floats, deg)
     alpha_t (n_eps*n_theta floats, row-major: outer eps, inner theta)
     alpha_n (n_eps*n_theta floats)
------------------------------------------------------------------------- */

void SurfCollideCLLMD::read_alpha_table(const char *fname)
{
  int me = comm->me;

  // read all non-comment tokens into a growable buffer

  double *buf = NULL;
  int nbuf = 0;

  if (me == 0) {
    FILE *fp = fopen(fname,"r");
    if (fp == NULL) {
      char str[256];
      sprintf(str,"Cannot open surf_collide cll/md alpha table file %s",fname);
      error->one(FLERR,str);
    }

    int maxbuf = 1024;
    buf = (double *) malloc(maxbuf*sizeof(double));
    char line[MAXLINE];
    while (fgets(line,MAXLINE,fp)) {
      char *s = line;
      while (*s == ' ' || *s == '\t') s++;
      if (*s == '#' || *s == '\0' || *s == '\n') continue;
      char *tok = strtok(line," \t\n\r");
      while (tok) {
        if (nbuf == maxbuf) {
          maxbuf *= 2;
          buf = (double *) realloc(buf,maxbuf*sizeof(double));
        }
        buf[nbuf++] = atof(tok);
        tok = strtok(NULL," \t\n\r");
      }
    }
    fclose(fp);

    if (nbuf < 2) error->one(FLERR,"surf_collide cll/md alpha table is empty");
    n_eps = (int) (buf[0] + 0.5);
    n_theta = (int) (buf[1] + 0.5);
    int expect = 2 + n_eps + n_theta + 2*n_eps*n_theta;
    if (n_eps < 2 || n_theta < 2 || nbuf < expect)
      error->one(FLERR,"surf_collide cll/md alpha table is malformed "
                 "(too few values for declared n_eps/n_theta)");
  }

  MPI_Bcast(&n_eps,1,MPI_INT,0,world);
  MPI_Bcast(&n_theta,1,MPI_INT,0,world);

  int ntab = n_eps*n_theta;
  eps_grid    = new double[n_eps];
  theta_grid  = new double[n_theta];
  alpha_t_tbl = new double[ntab];
  alpha_n_tbl = new double[ntab];

  if (me == 0) {
    int k = 2;
    for (int i = 0; i < n_eps; i++)   eps_grid[i]   = buf[k++];
    for (int j = 0; j < n_theta; j++) theta_grid[j] = buf[k++];
    for (int i = 0; i < ntab; i++)    alpha_t_tbl[i] = buf[k++];
    for (int i = 0; i < ntab; i++)    alpha_n_tbl[i] = buf[k++];
    free(buf);
  }

  MPI_Bcast(eps_grid,n_eps,MPI_DOUBLE,0,world);
  MPI_Bcast(theta_grid,n_theta,MPI_DOUBLE,0,world);
  MPI_Bcast(alpha_t_tbl,ntab,MPI_DOUBLE,0,world);
  MPI_Bcast(alpha_n_tbl,ntab,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   bilinear interpolation of alpha_t, alpha_n at (eps_eV, theta_deg).
   queries are clamped to the grid, so any extrapolation beyond the MD
   points must already be built into the table.
------------------------------------------------------------------------- */

void SurfCollideCLLMD::interpolate_alpha(double eps, double theta,
                                         double &a_t, double &a_n)
{
  if (eps < eps_grid[0]) eps = eps_grid[0];
  if (eps > eps_grid[n_eps-1]) eps = eps_grid[n_eps-1];
  if (theta < theta_grid[0]) theta = theta_grid[0];
  if (theta > theta_grid[n_theta-1]) theta = theta_grid[n_theta-1];

  int ie = 0;
  while (ie < n_eps-2 && eps > eps_grid[ie+1]) ie++;
  int it = 0;
  while (it < n_theta-2 && theta > theta_grid[it+1]) it++;

  double de = eps_grid[ie+1] - eps_grid[ie];
  double dt = theta_grid[it+1] - theta_grid[it];
  double fe = (de != 0.0) ? (eps - eps_grid[ie]) / de : 0.0;
  double ft = (dt != 0.0) ? (theta - theta_grid[it]) / dt : 0.0;
  if (fe < 0.0) fe = 0.0; if (fe > 1.0) fe = 1.0;
  if (ft < 0.0) ft = 0.0; if (ft > 1.0) ft = 1.0;

  int i00 = ie*n_theta + it;
  int i10 = (ie+1)*n_theta + it;
  int i01 = ie*n_theta + (it+1);
  int i11 = (ie+1)*n_theta + (it+1);

  double w00 = (1.0-fe)*(1.0-ft);
  double w10 = fe*(1.0-ft);
  double w01 = (1.0-fe)*ft;
  double w11 = fe*ft;

  a_t = w00*alpha_t_tbl[i00] + w10*alpha_t_tbl[i10] +
        w01*alpha_t_tbl[i01] + w11*alpha_t_tbl[i11];
  a_n = w00*alpha_n_tbl[i00] + w10*alpha_n_tbl[i10] +
        w01*alpha_n_tbl[i01] + w11*alpha_n_tbl[i11];

  // CLL needs acc in [0,1]; the table is clipped already, guard anyway
  if (a_t < 0.0) a_t = 0.0; if (a_t > 1.0) a_t = 1.0;
  if (a_n < 0.0) a_n = 0.0; if (a_n > 1.0) a_n = 1.0;
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   control flow as in stock cll::collide, reflection done by cll_md()
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideCLLMD::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;

  Particle::OnePart iorig;
  Particle::OnePart *jp = NULL;
  reaction = 0;
  int velreset = 0;

  if (isr >= 0) {
    if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    reaction = surf->sr[isr]->react(ip,isurf,norm,jp,velreset);
    if (reaction) surf->nreact_one++;
  }

  if (persurf_temperature) {
    tsurf = t_persurf[isurf];
    if (tsurf <= 0.0) error->one(FLERR,"Surf_collide tsurf <= 0.0");
  }

  if (ip) {
    if (!velreset) cll_md(ip,norm);
    if (modify->n_update_custom) {
      int i = ip - particle->particles;
      modify->update_custom(i,tsurf,tsurf,tsurf,vstream);
    }
  }
  if (jp) {
    if (!velreset) cll_md(jp,norm);
    if (modify->n_update_custom) {
      int j = jp - particle->particles;
      modify->update_custom(j,tsurf,tsurf,tsurf,vstream);
    }
  }

  if (reaction && modify->n_surf_react) {
    int i = -1;
    if (ip) i = ip - particle->particles;
    int j = -1;
    if (jp) j = jp - particle->particles;
    modify->surf_react(&iorig,i,j);
    if (jp && j < 0) {
      jp = NULL;
      particle->nlocal--;
    }
  }

  return jp;
}

/* ----------------------------------------------------------------------
   cll/md reflection: acc_n, acc_t from the MD table at (eps,theta)
------------------------------------------------------------------------- */

void SurfCollideCLLMD::cll_md(Particle::OnePart *p, double *norm)
{
  double tangent1[3],tangent2[3];
  Particle::Species *species = particle->species;
  int ispecies = p->ispecies;

  double *v = p->v;
  double dot = MathExtra::dot3(v,norm);
  double vmagsq = MathExtra::lensq3(v);
  double vrm, vperp, vtan1, vtan2;

  // incident energy (eV) and angle to the local normal (deg)

  double mass = species[ispecies].mass;
  double eps_eV = 0.5 * mass * vmagsq / EV2J;
  double theta_deg = 0.0;
  if (vmagsq > 0.0) {
    double costheta = fabs(dot) / sqrt(vmagsq);
    if (costheta > 1.0) costheta = 1.0;
    theta_deg = acos(costheta) * 180.0 / MY_PI;
  }

  double acc_n, acc_t;
  interpolate_alpha(eps_eV, theta_deg, acc_t, acc_n);

  // tangential basis

  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  if (MathExtra::lensq3(tangent1) == 0.0) {
    tangent2[0] = random->uniform();
    tangent2[1] = random->uniform();
    tangent2[2] = random->uniform();
    MathExtra::cross3(norm,tangent2,tangent1);
  }

  MathExtra::norm3(tangent1);
  MathExtra::cross3(norm,tangent1,tangent2);

  double tan1 = MathExtra::dot3(v,tangent1);

  vrm = sqrt(2.0*update->boltz * tsurf / mass);

  // normal component, acc_n from the table

  double r_1 = sqrt(-acc_n*log(random->uniform()));
  double theta_1 = MY_2PI * random->uniform();
  double dot_norm = dot/vrm * sqrt(1-acc_n);
  vperp = vrm * sqrt(r_1*r_1 + dot_norm*dot_norm + 2*r_1*dot_norm*cos(theta_1));

  // tangential components, acc_t from the table

  double r_2 = sqrt(-acc_t*log(random->uniform()));
  double theta_2 = MY_2PI * random->uniform();
  double vtangent = tan1/vrm * sqrt(1-acc_t);
  vtan1 = vrm * (vtangent + r_2*cos(theta_2));
  vtan2 = vrm * r_2 * sin(theta_2);

  v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
  v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
  v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

  // rotational energy, constant acc_rot

  if (!sparta->collide || sparta->collide->rotstyle == NONE ||
      species[ispecies].rotdof < 2) p->erot = 0.0;

  else {
    double erot_mag = sqrt(p->erot*(1-acc_rot)/(update->boltz*tsurf));
    double r_rot,cos_theta_rot,A_rot,X_rot;
    if (species[ispecies].rotdof == 2) {
      r_rot = sqrt(-acc_rot*log(random->uniform()));
      cos_theta_rot = cos(MY_2PI*random->uniform());
    }
    else {
      A_rot = 0;
      while (A_rot < random->uniform()) {
        X_rot = 4*random->uniform();
        A_rot = 2.71828182845904523536028747*X_rot*X_rot*exp(-X_rot*X_rot);
      }
      r_rot = sqrt(acc_rot)*X_rot;
      cos_theta_rot = 2*random->uniform() - 1;
    }
    p->erot = update->boltz * tsurf *
      (r_rot*r_rot + erot_mag*erot_mag + 2*r_rot*erot_mag*cos_theta_rot);
  }

  // vibrational energy, constant acc_vib

  int vibdof = species[ispecies].vibdof;
  double r_vib, cos_theta_vib, A_vib, X_vib, evib_mag;

  if (!sparta->collide || sparta->collide->vibstyle == NONE || vibdof < 2)
    p->evib = 0.0;

  else if (sparta->collide->vibstyle == DISCRETE && vibdof == 2) {
    double evib_star =
      -log(1 - random->uniform() *
           (1 - exp(-update->boltz*species[ispecies].vibtemp[0])));
    double evib_val = p->evib + evib_star;
    evib_mag = sqrt(evib_val*(1-acc_vib)/(update->boltz*tsurf));
    r_vib = sqrt(-acc_vib*log(random->uniform()));
    cos_theta_vib = cos(MY_2PI*random->uniform());
    evib_val = update->boltz * tsurf *
      (r_vib*r_vib + evib_mag*evib_mag + 2*r_vib*evib_mag*cos_theta_vib);
    int ivib =  evib_val / (update->boltz*species[ispecies].vibtemp[0]);
    p->evib = ivib * update->boltz * species[ispecies].vibtemp[0];
  }

  else if (sparta->collide->vibstyle == SMOOTH || vibdof >= 2) {
    evib_mag = sqrt(p->evib*(1-acc_vib)/(update->boltz*tsurf));
    if (vibdof == 2) {
      r_vib = sqrt(-acc_vib*log(random->uniform()));
      cos_theta_vib = cos(MY_2PI*random->uniform());
    } else {
      A_vib = 0;
      while (A_vib < random->uniform()) {
        X_vib = 4*random->uniform();
        A_vib = 2.71828182845904523536028747*X_vib*X_vib*exp(-X_vib*X_vib);
      }
      r_vib = sqrt(acc_vib)*X_vib;
      cos_theta_vib = 2*random->uniform() - 1;
    }
    p->evib = update->boltz * tsurf *
      (r_vib*r_vib + evib_mag*evib_mag + 2*r_vib*evib_mag*cos_theta_vib);
  }
}
