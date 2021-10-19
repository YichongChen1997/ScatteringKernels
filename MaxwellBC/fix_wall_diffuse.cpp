/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <cstring>
#include <cmath>            //SD
#include "fix_wall_diffuse.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "force.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

FixWallDiffuse::FixWallDiffuse(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nwall(0)
{
  if (narg < 6) error->all(FLERR,"Illegal fix wall/diffuse command");

  dynamic_group_allow = 1;

  // parse args

  nwall = 0;
  int scaleflag = 1;

  temperature = utils::numeric(FLERR,arg[3],false,lmp);
  alpha = utils::numeric(FLERR,arg[4],false,lmp);
  molMass = utils::numeric(FLERR,arg[5],false,lmp);
  int iarg = 6;
  while (iarg < narg) {
    if ((strcmp(arg[iarg],"xlo") == 0) || (strcmp(arg[iarg],"xhi") == 0) ||
        (strcmp(arg[iarg],"ylo") == 0) || (strcmp(arg[iarg],"yhi") == 0) ||
        (strcmp(arg[iarg],"zlo") == 0) || (strcmp(arg[iarg],"zhi") == 0)) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/diffuse command");

      int newwall;
      if (strcmp(arg[iarg],"xlo") == 0) newwall = XLO;
      else if (strcmp(arg[iarg],"xhi") == 0) newwall = XHI;
      else if (strcmp(arg[iarg],"ylo") == 0) newwall = YLO;
      else if (strcmp(arg[iarg],"yhi") == 0) newwall = YHI;
      else if (strcmp(arg[iarg],"zlo") == 0) newwall = ZLO;
      else if (strcmp(arg[iarg],"zhi") == 0) newwall = ZHI;

      for (int m = 0; (m < nwall) && (m < 6); m++)
        if (newwall == wallwhich[m])
          error->all(FLERR,"Wall defined twice in fix wall/diffuse command");

      wallwhich[nwall] = newwall;
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
        wallstyle[nwall] = EDGE;
        int dim = wallwhich[nwall] / 2;
        int side = wallwhich[nwall] % 2;
        if (side == 0) coord0[nwall] = domain->boxlo[dim];
        else coord0[nwall] = domain->boxhi[dim];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        wallstyle[nwall] = VARIABLE;
        int n = strlen(&arg[iarg+1][2]) + 1;
        varstr[nwall] = new char[n];
        strcpy(varstr[nwall],&arg[iarg+1][2]);
      } else {
        wallstyle[nwall] = CONSTANT;
        coord0[nwall] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      }

      nwall++;
      iarg += 2;

    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal wall/diffuse command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix wall/diffuse command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall/diffuse command");
  }

  // error check
  if (nwall == 0) error->all(FLERR,"Illegal fix wall command");

  for (int m = 0; m < nwall; m++) {
    if ((wallwhich[m] == XLO || wallwhich[m] == XHI) && domain->xperiodic)
      error->all(FLERR,"Cannot use fix wall/diffuse in periodic dimension");
    if ((wallwhich[m] == YLO || wallwhich[m] == YHI) && domain->yperiodic)
      error->all(FLERR,"Cannot use fix wall/diffuse in periodic dimension");
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->zperiodic)
      error->all(FLERR,"Cannot use fix wall/diffuse in periodic dimension");
  }

  for (int m = 0; m < nwall; m++)
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->dimension == 2)
      error->all(FLERR,
                 "Cannot use fix wall/diffuse zlo/zhi for a 2d simulation");

  // scale factors for CONSTANT and VARIABLE walls
  int flag = 0;
  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] != EDGE) flag = 1;

  if (flag) {
    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    }
    else xscale = yscale = zscale = 1.0;

    for (int m = 0; m < nwall; m++) {
      if (wallstyle[m] != CONSTANT) continue;
      if (wallwhich[m] < YLO) coord0[m] *= xscale;
      else if (wallwhich[m] < ZLO) coord0[m] *= yscale;
      else coord0[m] *= zscale;
    }
  }

  // set varflag if any wall positions are variable

  varflag = 0;
  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] == VARIABLE) varflag = 1;
}

/* ---------------------------------------------------------------------- */

FixWallDiffuse::~FixWallDiffuse()
{
  if (copymode) return;

  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] == VARIABLE) delete [] varstr[m];
}

/* ---------------------------------------------------------------------- */

int FixWallDiffuse::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallDiffuse::init()
{
  for (int m = 0; m < nwall; m++) {
    if (wallstyle[m] != VARIABLE) continue;
    varindex[m] = input->variable->find(varstr[m]);
    if (varindex[m] < 0)
      error->all(FLERR,"Variable name for fix wall/diffuse does not exist");
    if (!input->variable->equalstyle(varindex[m]))
      error->all(FLERR,"Variable for fix wall/diffuse is invalid style");
  }

  int nrigid = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;

  if (nrigid && comm->me == 0)
    error->warning(FLERR,"Should not allow rigid bodies to bounce off "
                   "diffusing walls");
}

/* ---------------------------------------------------------------------- */

void FixWallDiffuse::post_integrate()
{
  int i,dim,side;
  double coord;
  double rnd, theta, rho, refVelocity = 1e-10/1e-15, vx,vy,vz; //**** MB
  //double massOne;          //SD
  int dimOne, dimTwo;      //SD
  double term = sqrt(1.38E-23/molMass*temperature)/refVelocity; //**** MB
  bool diffuse = false;
  
  // coord = current position of wall
  // evaluate variable if necessary, wrap with clear/add

  double **x = atom->x;
  double **v = atom->v;
  //double *mass = atom->mass;        //SD
  //double *rmass = atom->rmass;      //SD
  //int *type = atom->type;           //SD
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  //double boltz = force->boltz;      //SD : Boltzmann constant in a unit that is specified in the lammps script.              
    
  if (varflag) modify->clearstep_compute();

  for (int m = 0; m < nwall; m++) {
    if (wallstyle[m] == VARIABLE) {
      coord = input->variable->compute_equal(varindex[m]);
      if (wallwhich[m] < YLO) coord *= xscale;
      else if (wallwhich[m] < ZLO) coord *= yscale;
      else coord *= zscale;
    } else coord = coord0[m];

    dim = wallwhich[m] / 2;
    side = wallwhich[m] % 2;
	
	if (dim == 0) {   // SD: to handle the velocity components other than the component normal to the surface. 
		dimOne = 1;
		dimTwo = 2;
	} else if (dim == 1) {
		dimOne = 0;
		dimTwo = 2;
	} else {
		dimOne = 0;
		dimTwo = 1;
	}		

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (side == 0) {
          if (x[i][dim] < coord) {
			  
              
             //**** MB start add *****
             
             //if (rmass) massOne = rmass[i];  
             // else massOne = mass[type[i]]; 
             /*SD: 'mass' is the per atom_type mass, while 'rmass' is the per atom_index mass. Most atom_style use per-type mass,
             but some atom_style (e.g. 'sphere') use per-atom rmass. */

            rnd =  0 + ( (double)rand() / RAND_MAX) * (1 - 0); // MB: random real number from 0 to 1
            diffuse = false;
            
            if (alpha == 1) diffuse = true;
            else if(rnd < alpha) diffuse = true;       
            
            if (diffuse)
            { 
                rnd =  0 + ( (double)rand() / RAND_MAX) * (1 - 0);  // v 
                theta = 2*3.141592*rnd;
                rnd = 0 + ( (double)rand() / RAND_MAX) * (1 - 0);   // u
                rho = sqrt( -2.0*log(rnd) );
                v[i][dimOne] = rho*cos(theta)*term;  // vx  *sqrt(boltz/massOne*temperature);
                v[i][dimTwo] = rho*sin(theta)*term;  // vy
                rnd =  0 + ( (double)rand() / RAND_MAX) * (1 - 0);  // j
                v[i][dim] = sqrt( -2.0*log(rnd) )*term;  //vz 
				
				        x[i][dim] = coord + (coord - x[i][dim]);
            }
            else // specular 
            {
                v[i][dim] = -v[i][dim];
			          x[i][dim] = coord + (coord - x[i][dim]);
            }

            //**** MB: end add --- now integrate vx,vy,vz into LAMMPS ****
            
            }
        } else {
          if (x[i][dim] > coord) {
			  
			//if (rmass) massOne = rmass[i];  
            //else massOne = mass[type[i]];   
                      
            rnd =  0 + ( (double)rand() / RAND_MAX) * (1 - 0); // MB: random real number from 0 to 1
			  
            diffuse = false;
            
            if (alpha == 1) diffuse = true;
            else if(rnd < alpha) diffuse = true;
            
            if (diffuse)
            {
                rnd =  0 + ( (double)rand() / RAND_MAX) * (1 - 0);  
                theta = 2*3.141592*rnd;
                rnd =  0 + ( (double)rand() / RAND_MAX) * (1 - 0); 
                rho = sqrt( -2.0*log(rnd) );
                v[i][dimOne] = rho*cos(theta)*term;       // vy  *sqrt(boltz/massOne*temperature);
                v[i][dimTwo] = rho*sin(theta)*term;       // vx
                rnd =  0 + ( (double)rand() / RAND_MAX) * (1 - 0);
                v[i][dim] = -sqrt( -2.0*log(rnd) )*term;  // vz
				
				        x[i][dim] = coord - (x[i][dim] - coord);
            }
            else
            {
                v[i][dim] = -v[i][dim];
			          x[i][dim] = coord - (x[i][dim] - coord);
            }                     
          }
        }
      }
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}
