/*!
  \file threed_singlefluid_fluxes.h
  \brief Riemann solvers for 3D single fluid Euler
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup fluxes
*/
#ifndef THREED_SINGLEFLUID_FLUXES_H
#define THREED_SINGLEFLUID_FLUXES_H
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
arch_device void threed_singlefluid_rusanov(scalar rhoL,
					    scalar rhoR,
					    scalar vxL,
					    scalar vxR,
					    scalar vyL,
					    scalar vyR,
					    scalar vzL,
					    scalar vzR,
					    scalar EtL,
					    scalar EtR,
					    scalar nx,
					    scalar ny,
					    scalar nz,
					    scalar* F, scalar* ncterm){
  
  scalar vnL = vxL*nx + vyL*ny + vzL*nz;
  scalar vnR = vxR*nx + vyR*ny + vzR*nz;
  scalar gamma = constants::GLOBAL_GAMMA;
  scalar pL = (gamma-1)*(EtL - 0.5*rhoL*(vxL*vxL + vyL*vyL + vzL*vzL));
  scalar pR = (gamma-1)*(EtR - 0.5*rhoR*(vxR*vxR + vyR*vyR + vzR*vzR));
  scalar aL = sqrt((gamma*pL)/rhoL);
  scalar aR = sqrt((gamma*pR)/rhoR);
  scalar HL = (EtL + pL)/rhoL;
  scalar HR = (EtR + pR)/rhoR;

  // Find the maximum eigenvalue
  scalar maxvap = MAX(fabs(vnL)+aL, fabs(vnR)+aR);

  //first: fx = rho*u; fy = rho*v; fz = rho*w
  F[0] = 0.5*(flux_ab(rhoL,vnL) + flux_ab(rhoR,vnR)
	      -maxvap*(rhoR-rhoL));

  //second:
  F[1] = 0.5*(flux_apb(rhoL*vnL*vxL , pL*nx)  + flux_apb(rhoR*vnR*vxR , pR*nx)
	      -maxvap*(rhoR*vxR - rhoL*vxL));

  //third:
  F[2] = 0.5*(flux_apb(rhoL*vnL*vyL , pL*ny)  + flux_apb(rhoR*vnR*vyR , pR*ny)
	      -maxvap*(rhoR*vyR - rhoL*vyL));

  //fourth:
  F[3] = 0.5*(flux_apb(rhoL*vnL*vzL , pL*nz)  + flux_apb(rhoR*vnR*vzR , pR*nz)
	      -maxvap*(rhoR*vzR - rhoL*vzL));

  //fifth: fx = rho*u*H; fy = rho*v*H; fz = rho*w*H
  F[4] = 0.5*(flux_abc(rhoL , vnL , HL)+flux_abc(rhoR , vnR , HR)
	      -maxvap*(EtR-EtL));

} // end Rusanov function
#endif //Rusanov endif

#ifdef CEN
arch_device void threed_singlefluid_central(scalar rhoL,
					    scalar rhoR,
					    scalar vxL,
					    scalar vxR,
					    scalar vyL,
					    scalar vyR,
					    scalar vzL,
					    scalar vzR,
					    scalar EtL,
					    scalar EtR,
					    scalar nx,
					    scalar ny,
					    scalar nz,
					    scalar* F, scalar* ncterm){

  //Just like Rusanov except we remove all maxvap terms
  scalar vnL = vxL*nx + vyL*ny + vzL*nz;
  scalar vnR = vxR*nx + vyR*ny + vzR*nz;
  scalar gamma = constants::GLOBAL_GAMMA;
  scalar pL = (gamma-1)*(EtL - 0.5*rhoL*(vxL*vxL + vyL*vyL + vzL*vzL));
  scalar pR = (gamma-1)*(EtR - 0.5*rhoR*(vxR*vxR + vyR*vyR + vzR*vzR));
  //scalar aL = sqrt((gamma*pL)/rhoL);
  //scalar aR = sqrt((gamma*pR)/rhoR);
  scalar HL = (EtL + pL)/rhoL;
  scalar HR = (EtR + pR)/rhoR;

  //first: fx = rho*u; fy = rho*v; fz = rho*w
  F[0] = 0.5*(flux_ab(rhoL,vnL) + flux_ab(rhoR,vnR));

  //second:
  F[1] = 0.5*(flux_apb(rhoL*vnL*vxL , pL*nx)  + flux_apb(rhoR*vnR*vxR , pR*nx));

  //third:
  F[2] = 0.5*(flux_apb(rhoL*vnL*vyL , pL*ny)  + flux_apb(rhoR*vnR*vyR , pR*ny));

  //fourth:
  F[3] = 0.5*(flux_apb(rhoL*vnL*vzL , pL*nz)  + flux_apb(rhoR*vnR*vzR , pR*nz));

  //fifth: fx = rho*u*H; fy = rho*v*H; fz = rho*w*H
  F[4] = 0.5*(flux_abc(rhoL , vnL , HL)+flux_abc(rhoR , vnR , HR));
	      
} // end central function
#endif //central flux endif

//*****************************************************************************
//* -- SLAU Flux Function ---
//*
//* Recently proposed member of the AUSM family. I'm using the SLAU 2
//* approach described in "Towards shock-stable and accurate hypersonic
//* heating computations: A new pressure flux for AUSM-family schemes"
//* by Kitamura and Shima
//* 
//* I was orignially interested in AUSM for low-Mach performance,
//* but SLAU2 approach is supposedely very good for all Mach numbers (nice bonus)
//*
//*****************************************************************************
#ifdef SLAU
arch_device void threed_singlefluid_slau2(scalar rhoL,
					  scalar rhoR,
					  scalar vxL,
					  scalar vxR,
					  scalar vyL,
					  scalar vyR,
					  scalar vzL,
					  scalar vzR,
					  scalar EtL,
					  scalar EtR,
					  scalar nx,
					  scalar ny,
					  scalar nz,
					  scalar* F, scalar* ncterm){
  
  scalar gamma = constants::GLOBAL_GAMMA;
  //Step 1: get pressure on each side
  scalar pL = (gamma-1)*(EtL - 0.5*rhoL*(vxL*vxL + vyL*vyL + vzL*vzL));
  scalar pR = (gamma-1)*(EtR - 0.5*rhoR*(vxR*vxR + vyR*vyR + vzR*vzR));

  //Step 2: Build left/right Psi vectors and the N vector
  scalar N[D+2];
  scalar PsiL[D+2];        scalar PsiR[D+2];
  PsiL[0] = 1.0;           PsiR[0] = 1.0;
  PsiL[1] = vxL;           PsiR[1] = vxR;
  PsiL[2] = vyL;           PsiR[2] = vyR;
  PsiL[3] = vzL;           PsiR[3] = vzR;
  PsiL[4] = (EtL+pL)/rhoL; PsiR[4] = (EtR+pR)/rhoR; //enthalpy
  N[0] = 0.0;
  N[1] = nx;
  N[2] = ny;
  N[3] = nz;
  N[4] = 0.0;

  //Step 3: Get left/right normal velocity
  scalar vnL = vxL*nx + vyL*ny + vzL*nz;
  scalar vnR = vxR*nx + vyR*ny + vzR*nz;

  //Step 4: Get speed of sound on each side,
  //then calculate arithmetic mean
  scalar cL = sqrt((gamma*pL)/rhoL);
  scalar cR = sqrt((gamma*pR)/rhoR);
  scalar cbar = 0.5*(cL+cR);

  //Step 5: Get left/right normal signed mach number
  scalar ML = vnL / cbar; scalar MR = vnR / cbar;

  //Step 6: Get the left/right f parameter
  scalar fL = 0.25*(ML+1.0) * (ML+1.0) * (2.0-ML);
  scalar fR = 0.25*(MR-1.0) * (MR-1.0) * (2.0+MR);
  if (fabs(ML)>1.0)
    {
      fL = 0.5*(1.0 + ML/fabs(ML)); //0.5 * (1+sign(ML))
    }
  if (fabs(MR)>1.0)
    {
      fR = 0.5*(1.0 - MR/fabs(MR)); //0.5 * (1-sign(MR))
    }

  //Step 7: Get the Mhat and Chi parameters
  scalar SS = sqrt(0.5*(vxL*vxL + vxR*vxR + vyL*vyL + vyR*vyR + vzL*vzL + vzR*vzR));
  scalar Mhat = fmin(1.0, SS/cbar);
  scalar Chi = (1.0-Mhat) * (1.0-Mhat);

  //Step 8: Get pressure flux(SLAU2 approach as opposed to SLAU)
  scalar rhobar = 0.5*(rhoL+rhoR);
  scalar pFlux = 0.5*(pL+pR) + 0.5*(fL-fR)*(pL-pR) + SS*(fL+fR-1.0)*rhobar*cbar;

  //Step 9: Get the g function
  scalar g = -1.0 * fmax(fmin(ML,0.0),-1.0) * fmin(fmax(MR,0.0),1.0);

  //Step 10: Some velocity stuff:
  scalar VnMag = (rhoL*fabs(vnL) + rhoR*fabs(vnR)) / (rhoL+rhoR);
  scalar VnMag_L = (1.0-g)*VnMag + g*fabs(vnL);
  scalar VnMag_R = (1.0-g)*VnMag + g*fabs(vnR);

  //Step 11: Get the momentum flux. 
  //Was initially unsure of sign in front of third term. Ran test with SOD
  //problem. Negative sign looks good, positive sign causes crash.
  //So, it will be (-Chi/cbar*(pR-pL))
  scalar mdot = 0.5*(rhoL*(vnL+VnMag_L) + rhoR*(vnR-VnMag_R) - Chi/cbar*(pR-pL));

  //Step 12: Assemble the flux
  scalar factorL = 0.5*(mdot+fabs(mdot)); scalar factorR = 0.5*(mdot-fabs(mdot));
  for (int j = 0; j < D+2; j++)
    {
      F[j] = flux_ab(factorL, PsiL[j]) + flux_ab(factorR, PsiR[j]) + pFlux*N[j];
    }
  //Done!
}
#endif //slau2 flux endif

#endif //threeD endif

#endif  //define header file endif
