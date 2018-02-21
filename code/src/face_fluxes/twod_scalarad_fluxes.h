/*!
  \file twod_singlefluid_fluxes.h
  \brief Riemann solvers for 2D single fluid Euler
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup fluxes
*/
#ifndef TWOD_SCALARAD_FLUXES_H
#define TWOD_SCALARAD_FLUXES_H
#ifdef TWOD
#include "scalar_def.h"
#include <math.h>
#include "macros.h"
#include "constants.h"
#include "basic_fluxes.h"
#include <stdio.h>

// Used to define dynamically mass fraction variables
#define YL(x) YL ##x
#define YR(x) YR ##x 

//***********************************************
/*
  Central flux function: just take average of the
  neighboring states. Be warned, this approach is generally
  frowned upon in the DG community
 */
//**********************************************
#ifdef CEN
arch_device void twod_scalarad_central(scalar rhoL,
					 scalar rhoR,
					 scalar nx,
					 scalar ny,
					 scalar* F, scalar* ncterm){

  scalar Vx = constants::GLOBAL_VX;
  scalar Vy = constants::GLOBAL_VY;
  scalar vnL = Vx*nx + Vy*ny;
  scalar vnR = Vx*nx + Vy*ny;

  //first: fx = rho*u; fy = rho*v
  F[0] = 0.5*(flux_ab(rhoL,vnL) + flux_ab(rhoR,vnR));
	      
} // end central function
#endif

//***********************************************
/*
  Upwind flux function: Perfect upwinding, analogous
  to the highly involved Godunov Riemann solution
  for the Euler equations.
 */
//**********************************************
#ifdef UPW
arch_device void twod_scalarad_upwind(scalar rhoL,
					 scalar rhoR,
					 scalar nx,
					 scalar ny,
					 scalar* F, scalar* ncterm){

  //assumption: nx,ny is out of L element and in to R element
  scalar Vx = constants::GLOBAL_VX;
  scalar Vy = constants::GLOBAL_VY;
  scalar vnL = Vx*nx + Vy*ny;
  scalar vnR = Vx*nx + Vy*ny;
  scalar maxvap = fabs(Vx);
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
  scalar Qx = Vx*0.5 * (rhoL*(1.0+sgn_E) + rhoR*(1.0-sgn_E)) * nx;
  scalar Qy = Vy*0.5 * (rhoL*(1.0+sgn_F) + rhoR*(1.0-sgn_F)) * ny;
  F[0] = Qx + Qy;
  

} // end central function
#endif

//*****************************************************************************
//* --- Rusanov's Flux Function ---
//*
//* V. V. Rusanov, Calculation of Interaction of Non-Steady Shock Waves with
//* Obstacles, J. Comput. Math. Phys. USSR, 1, pp. 267-279, 1961.
//*
//*****************************************************************************
#ifdef RUS
arch_device void twod_scalarad_rusanov(scalar rhoL,
					 scalar rhoR,
					 scalar nx,
					 scalar ny,
					 scalar* F, scalar* ncterm){

  scalar Vx = constants::GLOBAL_VX;
  scalar Vy = constants::GLOBAL_VY;
  scalar vnL = Vx*nx + Vy*ny;
  scalar vnR = Vx*nx + Vy*ny;

  // Find the maximum eigenvalue
  scalar maxvap = fabs(vnL);

  //first: fx = rho*u; fy = rho*v
  F[0] = 0.5*(flux_ab(rhoL,vnL) + flux_ab(rhoR,vnR)
	      -maxvap*(rhoR-rhoL));

} // end Rusanov function

#endif
#endif
#endif
