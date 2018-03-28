/*!
  \file twod_singlefluid_fluxes.h
  \brief Riemann solvers for 2D single fluid Euler
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup fluxes
*/
#ifndef THREED_SCALARAD_FLUXES_H
#define THREED_SCALARAD_FLUXES_H
#ifdef THREED
#include "scalar_def.h"
#include <math.h>
#include "macros.h"
#include "constants.h"
#include "basic_fluxes.h"
#include <stdio.h>

// Used to define dynamically mass fraction variables
#define YL(x) YL ##x
#define YR(x) YR ##x 

//*****************************************************************************
//* --- Rusanov's Flux Function ---
//*
//* V. V. Rusanov, Calculation of Interaction of Non-Steady Shock Waves with
//* Obstacles, J. Comput. Math. Phys. USSR, 1, pp. 267-279, 1961.
//*
//*****************************************************************************
#ifdef RUS
arch_device void threed_scalarad_rusanov(scalar rhoL,
					 scalar rhoR,
					 scalar nx,
					 scalar ny,
					 scalar nz,
					 scalar* F, scalar* ncterm){
  
  scalar Vx = constants::GLOBAL_VX;
  scalar Vy = constants::GLOBAL_VY;
  scalar Vz = constants::GLOBAL_VZ;
  
  scalar vnL = Vx*nx + Vy*ny + Vz*nz;
  scalar vnR = Vx*nx + Vy*ny + Vz*nz;
  
  // Find the maximum eigenvalue
  scalar maxvap = fabs(vnL);

  //first: fx = rho*u; fy = rho*v
  F[0] = 0.5*(flux_ab(rhoL,vnL) + flux_ab(rhoR,vnR)
	      -maxvap*(rhoR-rhoL));

} // end Rusanov function
#endif //Rusanov endif
#ifdef UPW
arch_device void threed_scalarad_upwind(scalar rhoL,
					scalar rhoR,
					scalar nx,
					scalar ny,
					scalar nz,
					scalar* F, scalar* ncterm){

  //assumption: nx,ny,nz is out of L element and in to R element
  scalar Vx = constants::GLOBAL_VX;
  scalar Vy = constants::GLOBAL_VY;
  scalar Vz = constants::GLOBAL_VZ;
  //scalar vnL = Vx*nx + Vy*ny + Vz*nz;
  //scalar vnR = Vx*nx + Vy*ny + Vz*nz;
  scalar sgn_E = 1.0;//sign(Vx*nx);
  if (Vx*nx < 0) 
    {
      sgn_E = -1.0;
    }
  scalar sgn_F = 1.0;
  if (Vy*ny < 0)
    {
      sgn_F = -1.0;
    }
  scalar sgn_G = 1.0;
  if (Vz*nz < 0)
    {
      sgn_G = -1.0;
    }
  scalar Qx = Vx*0.5 * (rhoL*(1.0+sgn_E) + rhoR*(1.0-sgn_E)) * nx;
  scalar Qy = Vy*0.5 * (rhoL*(1.0+sgn_F) + rhoR*(1.0-sgn_F)) * ny;
  scalar Qz = Vz*0.5 * (rhoL*(1.0+sgn_G) + rhoR*(1.0-sgn_G)) * nz;
  F[0] = Qx + Qy + Qz;
  

} // end upwind function
#endif //upwind flux endif
#ifdef CEN
arch_device void threed_scalarad_central(scalar rhoL,
					 scalar rhoR,
					 scalar nx,
					 scalar ny,
					 scalar nz,
					 scalar* F, scalar* ncterm){
  
  scalar Vx = constants::GLOBAL_VX;
  scalar Vy = constants::GLOBAL_VY;
  scalar Vz = constants::GLOBAL_VZ;
  scalar vnL = Vx*nx + Vy*ny + Vz*nz;
  scalar vnR = Vx*nx + Vy*ny + Vz*nz;

  //first: fx = rho*u; fy = rho*v
  F[0] = 0.5*(flux_ab(rhoL,vnL) + flux_ab(rhoR,vnR));
	      
} // end central function
#endif //central flux endif

#endif //threeD endif

#endif  //define header file endif
