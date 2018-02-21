/*!
  \file oned_singlefluid_fluxes.h
  \brief Riemann solvers for 1D Euler equations (singe fluids)
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup fluxes
*/
#ifndef ONED_SINGLEFLUID_FLUXES_H
#define ONED_SINGLEFLUID_FLUXES_H
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
arch_device void oned_singlefluid_rusanov(scalar rhoL,
					  scalar rhoR,
					  scalar vxL,
					  scalar vxR,
					  scalar EtL,
					  scalar EtR,
					  scalar nx,
					  scalar* F, scalar* ncterm){

  scalar gamma = constants::GLOBAL_GAMMA;
  scalar pL = (gamma-1)*(EtL - 0.5*rhoL*vxL*vxL);
  scalar pR = (gamma-1)*(EtR - 0.5*rhoR*vxR*vxR);

  scalar aL = sqrt((gamma*pL)/rhoL);
  scalar aR = sqrt((gamma*pR)/rhoR);
 
  // Find the maximum eigenvalue
  scalar maxvap = MAX(fabs(vxL)+aL, fabs(vxR)+aR);


  //first: fx = rho*u; 
  F[0] = 0.5*((flux_ab(rhoL,vxL) + flux_ab(rhoR,vxR))*nx
	      -maxvap*(rhoR-rhoL));

  //second: fx = rho*u*u+Bx*Bx+Pbar; 
  F[1] = 0.5*((flux_ab2pc(rhoL,vxL,pL)  + flux_ab2pc(rhoR,vxR,pR))*nx
	      -maxvap*(rhoR*vxR-rhoL*vxL));

  //third: fx = EtplusP*u; 
  F[2] = 0.5*((flux_ab(EtL+pL,vxL) + flux_ab(EtR+pR,vxR))*nx
	      -maxvap*(EtR-EtL));

} // end Rusanov function


//*****************************************************************************
//* --- HLL Flux Function ---
//*
//* A. Harten, P. D. Lax, and B. van Leer,On Upstream Differencing and
//* Godunov-Type Schemes for Hyperbolic Conservation Laws, SIAM Review,
//* 25(1), pp. 35-61, 1983.
//*
//* With default wave speeds evaluated by Davis's method:
//* Toro 10.48
//*
//* Also available wave speeds evaluated by Davis and Einfeldt's method:
//* Toro 10.49
//*
//* Also available wave speeds evaluated by Einfeldt's method:
//* Toro 10.52
//* B. Einfeldt, On Godunov-Type Methods for Gas Dynamics, SIAM Journal of
//* Numerical Analysis, 25(2), pp. 294-318, 1988.
//*
//*****************************************************************************
#elif HLL
arch_device void oned_singlefluid_hll(scalar rhoL,
				     scalar rhoR,
				     scalar vxL,
				     scalar vxR,
				     scalar EtL,
				     scalar EtR,
				     scalar nx,
				     scalar* F, scalar* ncterm){

#ifdef USE_CPU
  printf("Not implemented... The code is going to crash.\n");
#endif

} // end HLL function


//*****************************************************************************
//* -- Roe's Flux Function ---
//*
//* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
//* Schemes, Journal of Computational Physics, 43, pp. 357-372.
//*
//*****************************************************************************
#elif ROE
arch_device void oned_singlefluid_roe(scalar rhoL,
				     scalar rhoR,
				     scalar vxL,
				     scalar vxR,
				     scalar EtL,
				     scalar EtR,
				     scalar nx,
				     scalar* F, scalar* ncterm){

  scalar gamma = constants::GLOBAL_GAMMA;
  scalar pL = (gamma-1)*(EtL - 0.5*rhoL*vxL*vxL);
  scalar pR = (gamma-1)*(EtR - 0.5*rhoR*vxR*vxR);
 
  scalar aL = sqrt((gamma*pL)/rhoL);
  scalar aR = sqrt((gamma*pR)/rhoR);
  
  // Compute Roe averages
  scalar RT    = sqrt(rhoR/rhoL);
  scalar rho   = RT*rhoL;
  scalar v     = (vxL+RT*vxR)/(1+RT);
  scalar HL    = (EtL + pL)/rhoL;
  scalar HR    = (EtR + pR)/rhoR;
  scalar H     = (HL+RT*HR)/(1+RT);
  scalar a     = sqrt((gamma-1)*(H-0.5*v*v));
  scalar Dp  = pR - pL;
  
  // Roe waves strengths
  scalar dV0 = (Dp - rho*a*(vxR-vxL))/(2*a*a);
  scalar dV1 = (rhoR-rhoL) - Dp/(a*a);
  scalar dV2 = (Dp + rho*a*(vxR-vxL))/(2*a*a);

  // Absolute value of Roe eigenvalues
  scalar ws0 = fabs(v-a);
  scalar ws1 = fabs(v);
  scalar ws2 = fabs(v+a);
  
  // Roe Right eigenvectors
  scalar R00 = 1;
  scalar R01 = v-a;
  scalar R02 = H-v*a;

  scalar R10 = 1;
  scalar R11 = v;
  scalar R12 = 0.5*v*v;

  scalar R20 = 1;
  scalar R21 = v+a;
  scalar R22 = H+v*a;

  // first: fx = rho*u
  F[0] = 0.5*(flux_ab(rhoL,vxL) + flux_ab(rhoR,vxR))*nx
    -0.5*(ws0*dV0*R00+
	  ws1*dV1*R10+
	  ws2*dV2*R20);

  //second: fx = rho*u*u+Bx*Bx+Pbar; 
  F[1] = 0.5*(flux_ab2pc(rhoL,vxL,pL)  + flux_ab2pc(rhoR,vxR,pR))*nx
    -0.5*(ws0*dV0*R01+
	  ws1*dV1*R11+
	  ws2*dV2*R21);

  //third: fx = EtplusP*u; 
  F[2] = 0.5*(flux_ab(EtL+pL,vxL) + flux_ab(EtR+pR,vxR))*nx
    -0.5*(ws0*dV0*R02+
	  ws1*dV1*R12+
	  ws2*dV2*R22);

} // end Roe function
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
#elif SLAU
arch_device void oned_singlefluid_slau2(scalar rhoL,
					  scalar rhoR,
					  scalar vxL,
					  scalar vxR,
					  scalar EtL,
					  scalar EtR,
					  scalar nx,
					  scalar* F, scalar* ncterm){

scalar gamma = constants::GLOBAL_GAMMA;
  //Step 1: get pressure on each side
  scalar pL = (gamma-1)*(EtL - 0.5*rhoL*(vxL*vxL));
  scalar pR = (gamma-1)*(EtR - 0.5*rhoR*(vxR*vxR));

  //Step 2: Build left/right Psi vectors and the N vector
  scalar N[D+2];
  scalar PsiL[D+2];        scalar PsiR[D+2];
  PsiL[0] = 1.0;           PsiR[0] = 1.0;
  PsiL[1] = vxL;           PsiR[1] = vxR;
  PsiL[2] = (EtL+pL)/rhoL; PsiR[2] = (EtR+pR)/rhoR; //enthalpy
  N[0] = 0.0;
  N[1] = nx;
  N[2] = 0.0;

  //Step 3: Get left/right normal velocity
  scalar vnL = vxL*nx;
  scalar vnR = vxR*nx;

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
  if (fabs(ML>1))
    {
      fL = 0.5*(1.0 + ML/fabs(ML)); //0.5 * (1+sign(ML))
    }
  if (fabs(MR>1))
    {
      fR = 0.5*(1.0 - MR/fabs(MR)); //0.5 * (1-sign(MR))
}

  //Step 7: Get the Mhat and Chi parameters
scalar SS = sqrt(0.5*(vxL*vxL + vxR*vxR));
  scalar Mhat = fmin(1.0, SS/cbar);
  scalar Chi = (1.0-Mhat) * (1.0-Mhat);

  //Step 8: Get pressure flux(SLAU2 approach as opposed to SLAU)
  scalar rhobar = 0.5*(rhoL+rhoR);
  scalar pFlux = 0.5*(pL+pR) + 0.5*(fL-fR)*(pL-pR) + SS*(fL+fR-1.0)*rhobar*cbar;

  //Step 9: Get the g function
  scalar g = -1.0 * fmax(fmin(ML,0.0),-1.0) * fmin(fmax(MR,0.0),1.0);

  //Step 10: Some velocity stuff:
  scalar VnMag = (rhoL*fabs(vnL) + rhoR*fabs(vnR)) / (rhoL+rhoR);
  scalar VnLMag = (1.0-g)*VnMag + g*fabs(vnL);
  scalar VnRMag = (1.0-g)*VnMag + g*fabs(vnR);

  //Step 11: Get the momentum flux. 
  //Was initially unsure of sign in front of third term. Ran test with SOD
  //problem. Negative sign looks good, positive sign causes crash.
  //So, it will be (-Chi/cbar*(pR-pL))
  scalar mdot = 0.5*(rhoL*(vnL+VnLMag) + rhoR*(vnR-VnRMag) - Chi/cbar*(pR-pL));

  //Step 12: Assemble the flux
  scalar factorL = 0.5*(mdot+fabs(mdot)); scalar factorR = 0.5*(mdot-fabs(mdot));
  for (int j = 0; j < D+2; j++)
    {
      F[j] = flux_ab(factorL, PsiL[j]) + flux_ab(factorR, PsiR[j]) + pFlux*N[j];
    }
  //Done!
} //end SLAU2 function
#endif
#endif
#endif
