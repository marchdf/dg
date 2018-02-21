/*!
  \file oned_singlefluid_fluxes.h
  \brief Riemann solvers for 1D Euler equations (singe fluids)
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup fluxes
*/
#ifndef ONED_SCALARAD_FLUXES_H
#define ONED_SCALARAD_FLUXES_H
#ifdef ONED
#include "scalar_def.h"
#include <math.h>
#include "macros.h"
#include "constants.h"
#include "basic_fluxes.h"
#include <stdio.h>

//*****************************************************************************
//* --- Rusanov's Flux Function ---
//*
//* V. V. Rusanov, Calculation of Interaction of Non-Steady Shock Waves with
//* Obstacles, J. Comput. Math. Phys. USSR, 1, pp. 267-279, 1961.
//*
//*****************************************************************************
#ifdef RUS
arch_device void oned_scalarad_rusanov(scalar rhoL,
					  scalar rhoR,
					  scalar nx,
					  scalar* F, scalar* ncterm){

  scalar Vx = constants::GLOBAL_VX;
  scalar maxvap = fabs(Vx);
  
  
  //first: fx = rho*Vx; 
  //  F[0] = 0.5*((flux_ab(rhoL,Vx) + flux_ab(rhoR,Vx))
  //	      -0.5*maxvap*fabs(rhoR-rhoL)) * nx;

  //Copied from singlefulid rusanov(Marc's work)
  //F[0] = 0.5*((flux_ab(rhoL,Vx) + flux_ab(rhoR,Vx))*nx
  //	      -maxvap*(rhoR-rhoL));

  //Modified from singlefulid rusanov(Marc's work)
  F[0] = 0.5*((flux_ab(rhoL,Vx) + flux_ab(rhoR,Vx))*nx
	      -2.0*maxvap*(rhoR-rhoL));

} // end Rusanov function
#endif

//Pure upwind flux: the analogue is the full Godunov solution in Euler equation case
#ifdef UPW

arch_device void oned_scalarad_upwind(scalar rhoL,
					  scalar rhoR,
					  scalar nx,
					  scalar* F, scalar* ncterm){

  //assumption: nx is out of L element and in to R element
  scalar Vx = constants::GLOBAL_VX;
  scalar maxvap = fabs(Vx);
  scalar sgn_E = 1.0;//sign(Vx*nx);
  if (Vx*nx < 0)
    {
      sgn_E = -1.0;
    }
  scalar Qx = Vx*0.5 * (rhoL*(1.0+sgn_E) + rhoR*(1.0-sgn_E)) * nx;
  F[0] = Qx;
  

} // end upwind flux
#endif

//Central flux: Take the average flux of the two sides
#ifdef CEN
arch_device void oned_scalarad_central(scalar rhoL,
					  scalar rhoR,
					  scalar nx,
					  scalar* F, scalar* ncterm){

  scalar Vx = constants::GLOBAL_VX;
  scalar maxvap = fabs(Vx);
  
  
  //first: fx = rho*Vx; 
  F[0] = 0.5*((flux_ab(rhoL,Vx) + flux_ab(rhoR,Vx)))*nx;

} // end Rusanov function
#endif //end central flux

#endif
#endif 
