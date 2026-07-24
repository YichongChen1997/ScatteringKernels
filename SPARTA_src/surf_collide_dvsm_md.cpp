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

   dvsm/md: non-separable surface collision by direct resampling of MD
   beam events. See surf_collide_dvsm_md.h for the model.
   Usage: surf_collide <id> dvsm/md Tsurf <event_library>
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "surf_collide_dvsm_md.h"
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

SurfCollideDVSMMD::SurfCollideDVSMMD(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  // surf_collide <id> dvsm/md Tsurf <event_library>

  if (narg != 4) error->all(FLERR,"Illegal surf_collide dvsm/md command");

  parse_tsurf(arg[2]);

  eps_grid = theta_grid = ptrap = events = NULL;
  nevents = offset = NULL;
  n_eps = n_theta = 0;
  read_event_library(arg[3]);

  vstream[0] = vstream[1] = vstream[2] = 0.0;

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideDVSMMD::~SurfCollideDVSMMD()
{
  delete random;
  delete [] eps_grid;
  delete [] theta_grid;
  delete [] nevents;
  delete [] ptrap;
  delete [] offset;
  delete [] events;
}

/* ---------------------------------------------------------------------- */

void SurfCollideDVSMMD::init()
{
  SurfCollide::init();
  check_tsurf();
}

/* ----------------------------------------------------------------------
   read the MD event library on rank 0, broadcast to all ranks.
   File format (# starts a comment, tokens may be split over lines):
     n_eps n_theta
     eps_grid (n_eps floats, eV)
     theta_grid (n_theta floats, deg)
     per bin (outer eps, inner theta):
       eps theta N_return p_trap
       N_return x (vx vy vz)   [m/s, beam frame: +z outward normal,
                                +x incident tangential]
------------------------------------------------------------------------- */

void SurfCollideDVSMMD::read_event_library(const char *fname)
{
  int me = comm->me;

  // read all non-comment tokens into a growable buffer

  double *buf = NULL;
  bigint nbuf = 0;
  bigint ntotal = 0;

  if (me == 0) {
    FILE *fp = fopen(fname,"r");
    if (fp == NULL) {
      char str[256];
      sprintf(str,"Cannot open surf_collide dvsm/md event library file %s",
              fname);
      error->one(FLERR,str);
    }

    bigint maxbuf = 65536;
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

    if (nbuf < 2) error->one(FLERR,"surf_collide dvsm/md event library "
                             "is empty");
    n_eps = (int) (buf[0] + 0.5);
    n_theta = (int) (buf[1] + 0.5);
    if (n_eps < 2 || n_theta < 2)
      error->one(FLERR,"surf_collide dvsm/md event library grid too small");
  }

  MPI_Bcast(&n_eps,1,MPI_INT,0,world);
  MPI_Bcast(&n_theta,1,MPI_INT,0,world);

  int nbin = n_eps*n_theta;
  eps_grid   = new double[n_eps];
  theta_grid = new double[n_theta];
  nevents    = new int[nbin];
  ptrap      = new double[nbin];
  offset     = new int[nbin];

  if (me == 0) {
    bigint k = 2;
    if (k + n_eps + n_theta > nbuf)
      error->one(FLERR,"surf_collide dvsm/md event library is malformed");
    for (int i = 0; i < n_eps; i++)   eps_grid[i]   = buf[k++];
    for (int j = 0; j < n_theta; j++) theta_grid[j] = buf[k++];

    // first pass: bin headers and counts, record where each bin starts

    bigint *evstart = new bigint[nbin];
    for (int ie = 0; ie < n_eps; ie++) {
      for (int it = 0; it < n_theta; it++) {
        int b = ie*n_theta + it;
        if (k + 4 > nbuf)
          error->one(FLERR,"surf_collide dvsm/md event library is malformed");
        double e_hdr = buf[k++];
        double t_hdr = buf[k++];
        int n = (int) (buf[k] + 0.5); k++;
        double pt = buf[k++];
        if (fabs(e_hdr - eps_grid[ie]) > 1.0e-6*eps_grid[ie] ||
            fabs(t_hdr - theta_grid[it]) > 1.0e-6*theta_grid[it])
          error->one(FLERR,"surf_collide dvsm/md event library bin header "
                     "does not match declared grid order");
        if (n < 0 || pt < 0.0 || pt > 1.0)
          error->one(FLERR,"surf_collide dvsm/md event library has invalid "
                     "bin count or p_trap");
        if (n == 0 && pt < 1.0)
          error->one(FLERR,"surf_collide dvsm/md event library bin has "
                     "no events but p_trap < 1");
        nevents[b] = n;
        ptrap[b] = pt;
        evstart[b] = k;
        k += 3*(bigint)n;
        if (k > nbuf)
          error->one(FLERR,"surf_collide dvsm/md event library is truncated");
      }
    }

    offset[0] = 0;
    for (int b = 1; b < nbin; b++) offset[b] = offset[b-1] + nevents[b-1];
    ntotal = offset[nbin-1] + nevents[nbin-1];

    events = new double[3*ntotal];
    for (int b = 0; b < nbin; b++) {
      bigint src = evstart[b];
      for (int i = 0; i < nevents[b]; i++) {
        double vx = buf[src++];
        double vy = buf[src++];
        double vz = buf[src++];
        if (vz <= 0.0)
          error->one(FLERR,"surf_collide dvsm/md event library has "
                     "non-outgoing event (vz <= 0)");
        events[3*(offset[b]+i)  ] = vx;
        events[3*(offset[b]+i)+1] = vy;
        events[3*(offset[b]+i)+2] = vz;
      }
    }
    delete [] evstart;
    free(buf);
  }

  MPI_Bcast(eps_grid,n_eps,MPI_DOUBLE,0,world);
  MPI_Bcast(theta_grid,n_theta,MPI_DOUBLE,0,world);
  MPI_Bcast(nevents,nbin,MPI_INT,0,world);
  MPI_Bcast(ptrap,nbin,MPI_DOUBLE,0,world);
  MPI_Bcast(offset,nbin,MPI_INT,0,world);
  MPI_Bcast(&ntotal,1,MPI_SPARTA_BIGINT,0,world);
  if (me != 0) events = new double[3*ntotal];
  MPI_Bcast(events,3*ntotal,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   control flow as in cll/md::collide, reflection done by dvsm()
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideDVSMMD::
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
    if (!velreset) dvsm(ip,norm);
    if (modify->n_update_custom) {
      int i = ip - particle->particles;
      modify->update_custom(i,tsurf,tsurf,tsurf,vstream);
    }
  }
  if (jp) {
    if (!velreset) dvsm(jp,norm);
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
   dvsm reflection: draw an outgoing velocity from the MD event library.
   Bin selection is stochastic bilinear, in ln(eps) for the geometric
   energy grid and linear in angle. Incident states outside the grid
   clamp to the edge bins, as the alpha table of cll/md does, so both
   boundary conditions have the same grid support. The outgoing speed is
   rescaled by sqrt(eps/eps_bin), which keeps the MD fractional energy
   exchange at off-grid incident energies.
------------------------------------------------------------------------- */

void SurfCollideDVSMMD::dvsm(Particle::OnePart *p, double *norm)
{
  double tangent1[3],tangent2[3];
  Particle::Species *species = particle->species;
  int ispecies = p->ispecies;

  double *v = p->v;
  double dot = MathExtra::dot3(v,norm);
  double vmagsq = MathExtra::lensq3(v);

  // tangential basis, tangent1 along the incident tangential component

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

  if (vmagsq == 0.0) {
    diffuse_emit(p,norm,tangent1,tangent2);
  } else {

    // incident energy (eV) and angle to the local normal (deg)

    double mass = species[ispecies].mass;
    double eps_eV = 0.5 * mass * vmagsq / EV2J;
    double costheta = fabs(dot) / sqrt(vmagsq);
    if (costheta > 1.0) costheta = 1.0;
    double theta_deg = acos(costheta) * 180.0 / MY_PI;

    // stochastic bilinear bin selection

    int ie;
    if (eps_eV <= eps_grid[0]) ie = 0;
    else if (eps_eV >= eps_grid[n_eps-1]) ie = n_eps-1;
    else {
      ie = 0;
      while (ie < n_eps-2 && eps_eV > eps_grid[ie+1]) ie++;
      double de = log(eps_grid[ie+1]) - log(eps_grid[ie]);
      double fe = (de != 0.0) ? (log(eps_eV) - log(eps_grid[ie])) / de : 0.0;
      if (random->uniform() < fe) ie++;
    }

    int it;
    if (theta_deg <= theta_grid[0]) it = 0;
    else if (theta_deg >= theta_grid[n_theta-1]) it = n_theta-1;
    else {
      it = 0;
      while (it < n_theta-2 && theta_deg > theta_grid[it+1]) it++;
      double dt = theta_grid[it+1] - theta_grid[it];
      double ft = (dt != 0.0) ? (theta_deg - theta_grid[it]) / dt : 0.0;
      if (random->uniform() < ft) it++;
    }

    int bin = ie*n_theta + it;

    if (nevents[bin] == 0 || random->uniform() < ptrap[bin]) {

      // trapping-desorption channel: thermal re-emission at Tsurf

      diffuse_emit(p,norm,tangent1,tangent2);

    } else {

      // impulsive channel: one MD event, rotated into the local frame

      int k = (int) (nevents[bin] * random->uniform());
      if (k >= nevents[bin]) k = nevents[bin]-1;
      const double *ev = &events[3*(offset[bin]+k)];

      double s = sqrt(eps_eV / eps_grid[ie]);
      double sy = (random->uniform() < 0.5) ? -1.0 : 1.0;
      double vt1 = s*ev[0];
      double vt2 = s*sy*ev[1];
      double vperp = s*ev[2];

      v[0] = vperp*norm[0] + vt1*tangent1[0] + vt2*tangent2[0];
      v[1] = vperp*norm[1] + vt1*tangent1[1] + vt2*tangent2[1];
      v[2] = vperp*norm[2] + vt1*tangent1[2] + vt2*tangent2[2];
    }
  }

  // the event library carries no internal state, so rot/vib relax fully
  // to Tsurf, the acc = 1 limit of the CLL expressions

  if (!sparta->collide || sparta->collide->rotstyle == NONE ||
      species[ispecies].rotdof < 2) p->erot = 0.0;
  else {
    if (species[ispecies].rotdof == 2)
      p->erot = -update->boltz * tsurf * log(random->uniform());
    else {
      double A_rot = 0.0, X_rot = 0.0;
      while (A_rot < random->uniform()) {
        X_rot = 4*random->uniform();
        A_rot = 2.71828182845904523536028747*X_rot*X_rot*exp(-X_rot*X_rot);
      }
      p->erot = update->boltz * tsurf * X_rot*X_rot;
    }
  }

  int vibdof = species[ispecies].vibdof;
  if (!sparta->collide || sparta->collide->vibstyle == NONE || vibdof < 2)
    p->evib = 0.0;
  else if (sparta->collide->vibstyle == DISCRETE && vibdof == 2) {
    double evib_val = -update->boltz * tsurf * log(random->uniform());
    int ivib = evib_val / (update->boltz*species[ispecies].vibtemp[0]);
    p->evib = ivib * update->boltz * species[ispecies].vibtemp[0];
  }
  else {
    if (vibdof == 2)
      p->evib = -update->boltz * tsurf * log(random->uniform());
    else {
      double A_vib = 0.0, X_vib = 0.0;
      while (A_vib < random->uniform()) {
        X_vib = 4*random->uniform();
        A_vib = 2.71828182845904523536028747*X_vib*X_vib*exp(-X_vib*X_vib);
      }
      p->evib = update->boltz * tsurf * X_vib*X_vib;
    }
  }
}

/* ----------------------------------------------------------------------
   thermal re-emission at Tsurf, velocity only.
   the caller handles internal energy
------------------------------------------------------------------------- */

void SurfCollideDVSMMD::diffuse_emit(Particle::OnePart *p, double *norm,
                                     double *tangent1, double *tangent2)
{
  Particle::Species *species = particle->species;
  int ispecies = p->ispecies;

  double vrm = sqrt(2.0*update->boltz * tsurf / species[ispecies].mass);
  double vperp = vrm * sqrt(-log(random->uniform()));

  double theta = MY_2PI * random->uniform();
  double vtangent = vrm * sqrt(-log(random->uniform()));
  double vtan1 = vtangent * sin(theta);
  double vtan2 = vtangent * cos(theta);

  double *v = p->v;
  v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
  v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
  v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];
}
