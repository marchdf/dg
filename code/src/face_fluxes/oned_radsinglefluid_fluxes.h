/*!
  \file oned_sradinglefluid_fluxes.h
  \brief Riemann solvers for 1D Euler equations ,singlefluid with Reisner artificial dissipation
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup fluxes
*/
#ifndef ONED_RADSINGLEFLUID_FLUXES_H
#define ONED_RADSINGLEFLUID_FLUXES_H
#ifdef ONED
#include "scalar_def.h"
#include <math.h>
#include "macros.h"
#include "constants.h"
#include "basic_fluxes.h"
#include <stdio.h>

#ifdef SLAU
arch_device void oned_radsinglefluid_slau2(scalar rhoL,
					   scalar rhoR,
					   scalar vxL,
					   scalar vxR,
					   scalar EtL,
					   scalar EtR,
					   scalar CL,
					   scalar CR,
					   scalar nx,
					   scalar* F, scalar* ncterm){
  //printf("SLAU2 Riemann solver: (rho,u,Et,C)_L  = (%f, %f, %f, %f), (rho,u,Et,C)_R  = (%f, %f, %f, %f)\n",rhoL,vxL,EtL,CL,rhoR,vxR,EtR,CR);
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
      //      printf("In Slau2 Riemann solver: j=%d, F=%f\n",j,F[j]);
    }
  //Step 12b for radsinglefluid: flux = 0 because I do all flux calculations in q_rad routine
  F[D+2] = 0.0;

  //Done!
} //end SLAU2 function
#endif
#endif
#endif
