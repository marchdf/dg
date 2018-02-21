/*!
  \file twod_singlefluid_fluxes.h
  \brief Riemann solvers for 2D single fluid Euler
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup fluxes
*/
#ifndef TWOD_RADSINGLEFLUID_FLUXES_H
#define TWOD_RADSINGLEFLUID_FLUXES_H
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

#ifdef ROE
arch_device void twod_radsinglefluid_roe(scalar rhoL,
				     scalar rhoR,
				     scalar vxL,
				     scalar vxR,
				     scalar vyL,
				     scalar vyR,
				     scalar EtL,
				     scalar EtR,
				     scalar nx,
				     scalar ny,
				     scalar* F, scalar* ncterm){

  scalar tx = -ny;
  scalar ty =  nx;
  scalar vnL   = vxL*nx+vyL*ny;
  scalar vnR   = vxR*nx+vyR*ny;
  scalar vtL   = vxL*tx+vyL*ty;
  scalar vtR   = vxR*tx+vyR*ty;
  scalar gamma = constants::GLOBAL_GAMMA;
  scalar pL = (gamma-1)*(EtL - 0.5*rhoL*(vxL*vxL+vyL*vyL));
  scalar pR = (gamma-1)*(EtR - 0.5*rhoR*(vxR*vxR+vyR*vyR));
  scalar aL = sqrt((gamma*pL)/rhoL);
  scalar aR = sqrt((gamma*pR)/rhoR);
  scalar HL = (EtL + pL)/rhoL;
  scalar HR = (EtR + pR)/rhoR;

  // Compute Roe averages
  scalar RT  = sqrt(rhoR/rhoL);
  scalar rho = RT*rhoL;
  scalar vx  = (vxL+RT*vxR)/(1+RT);
  scalar vy  = (vyL+RT*vyR)/(1+RT);
  scalar vn  = vx*nx+vy*ny; 
  scalar vt  = vx*tx+vy*ty; 
  scalar H   = (HL+RT*HR)/(1+RT);
  scalar a   = sqrt((gamma-1)*(H-0.5*(vx*vx+vy*vy)));
 
  // Roe waves strengths
  scalar drho = rhoR - rhoL;
  scalar dp   = pR - pL;
  scalar dvn =  vnR - vnL;
  scalar dvt =  vtR - vtL;
  scalar dV0 = (dp - rho*a*dvn )/(2*a*a);
  scalar dV1 = rho*dvt/a;
  scalar dV2 = drho - dp/(a*a);
  scalar dV3 = (dp + rho*a*dvn )/(2*a*a);

  // Absolute value of Roe eigenvalues
  scalar ws0 = fabs(vn-a);
  scalar ws1 = fabs(vn);
  scalar ws2 = fabs(vn);
  scalar ws3 = fabs(vn+a);


  // Roe Right eigenvectors
  scalar R00 = 1;
  scalar R01 = vx - a*nx;
  scalar R02 = vy - a*ny;
  scalar R03 = H - vn*a;

  scalar R10 = 0;
  scalar R11 = a*tx;
  scalar R12 = a*ty;
  scalar R13 = vt*a;

  scalar R20 = 1;
  scalar R21 = vx;
  scalar R22 = vy;
  scalar R23 = 0.5*(vx*vx+vy*vy);

  scalar R30 = 1;
  scalar R31 = vx + a*nx;
  scalar R32 = vy + a*ny;
  scalar R33 = H + vn*a;

  //first: fx = rho*u; fy = rho*v
  F[0] = 0.5*(flux_ab(rhoL,vnL) + flux_ab(rhoR,vnR))
    -0.5*(ws0*dV0*R00+
	  ws1*dV1*R10+
	  ws2*dV2*R20+
	  ws3*dV3*R30);

  //second: fx = rho*u*u+p; fy = rho*u*v
  F[1] = 0.5*(flux_apb(rhoL*vnL*vxL,pL*nx)  + flux_apb(rhoR*vnR*vxR,pR*nx))
    -0.5*(ws0*dV0*R01+
	  ws1*dV1*R11+
	  ws2*dV2*R21+
	  ws3*dV3*R31);

  //third: fx = rho*u*v; fy = rho*v*v+p
  F[2] = 0.5*(flux_apb(rhoL*vnL*vyL,pL*ny)  + flux_apb(rhoR*vnR*vyR,pR*ny))
    -0.5*(ws0*dV0*R02+
	  ws1*dV1*R12+
	  ws2*dV2*R22+
	  ws3*dV3*R32);
 
  //fourth: fx = rho*u*H; fy = rho*v*H;
  F[3] = 0.5*(flux_abc(rhoL,vnL,HL)+flux_abc(rhoR,vnR,HR))
    -0.5*(ws0*dV0*R03+
	  ws1*dV1*R13+
	  ws2*dV2*R23+
	  ws3*dV3*R33);

  F[4] = 0.0;
  F[5] = 0.0;
  
} // end Roe function
#endif

#ifdef RUS
arch_device void twod_radsinglefluid_rusanov(scalar rhoL,
					 scalar rhoR,
					 scalar vxL,
					 scalar vxR,
					 scalar vyL,
					 scalar vyR,
					 scalar EtL,
					 scalar EtR,
					 scalar nx,
					 scalar ny,
					 scalar* F, scalar* ncterm){

  scalar vnL = vxL*nx+vyL*ny;
  scalar vnR = vxR*nx+vyR*ny;
  scalar gamma = constants::GLOBAL_GAMMA;
  scalar pL = (gamma-1)*(EtL - 0.5*rhoL*(vxL*vxL+vyL*vyL));
  scalar pR = (gamma-1)*(EtR - 0.5*rhoR*(vxR*vxR+vyR*vyR));
  scalar aL = sqrt((gamma*pL)/rhoL);
  scalar aR = sqrt((gamma*pR)/rhoR);
  scalar HL = (EtL + pL)/rhoL;
  scalar HR = (EtR + pR)/rhoR;

  // Find the maximum eigenvalue
  scalar maxvap = MAX(fabs(vnL)+aL, fabs(vnR)+aR);

  //first: fx = rho*u; fy = rho*v
  F[0] = 0.5*(flux_ab(rhoL,vnL) + flux_ab(rhoR,vnR)
	      -maxvap*(rhoR-rhoL));

  //second: fx = rho*u*u+p; fy = rho*u*v
  F[1] = 0.5*(flux_apb(rhoL*vnL*vxL,pL*nx)  + flux_apb(rhoR*vnR*vxR,pR*nx)
	      -maxvap*(rhoR*vxR-rhoL*vxL));

  //third: fx = rho*u*v; fy = rho*v*v+p
  F[2] = 0.5*(flux_apb(rhoL*vnL*vyL,pL*ny)  + flux_apb(rhoR*vnR*vyR,pR*ny)
	      -maxvap*(rhoR*vyR-rhoL*vyL));

  //fourth: fx = rho*u*H; fy = rho*v*H;
  F[3] = 0.5*(flux_abc(rhoL,vnL,HL)+flux_abc(rhoR,vnR,HR)
	      -maxvap*(EtR-EtL));

  //The C-equation fluxes: Set to zero
  F[4] = 0.0;
  F[5] = 0.0;

} // end Rusanov function
#endif

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
arch_device void twod_radsinglefluid_slau2(scalar rhoL,
					   scalar rhoR,
					   scalar vxL,
					   scalar vxR,
					   scalar vyL,
					   scalar vyR,
					   scalar EtL,
					   scalar EtR,
					   scalar nx,
					   scalar ny,
					   scalar* F, scalar* ncterm){
  
  scalar gamma = constants::GLOBAL_GAMMA;
  //Step 1: get pressure on each side
  scalar pL = (gamma-1)*(EtL - 0.5*rhoL*(vxL*vxL + vyL*vyL));
  scalar pR = (gamma-1)*(EtR - 0.5*rhoR*(vxR*vxR + vyR*vyR));

  //Step 2: Build left/right Psi vectors and the N vector
  scalar N[D+2];
  scalar PsiL[D+2];        scalar PsiR[D+2];
  PsiL[0] = 1.0;           PsiR[0] = 1.0;
  PsiL[1] = vxL;           PsiR[1] = vxR;
  PsiL[2] = vyL;           PsiR[2] = vyR;
  PsiL[3] = (EtL+pL)/rhoL; PsiR[3] = (EtR+pR)/rhoR; //enthalpy
  N[0] = 0.0;
  N[1] = nx;
  N[2] = ny;
  N[3] = 0.0;

  //Step 3: Get left/right normal velocity
  scalar vnL = vxL*nx + vyL*ny;
  scalar vnR = vxR*nx + vyR*ny;

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
  scalar SS = sqrt(0.5*(vxL*vxL + vxR*vxR + vyL*vyL + vyR*vyR));
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

  //Just one change from usual singlefluid case: set face fluxes for C variable to zero
  F[D+2] = 0.0;
  F[D+3] = 0.0;
  //Done!
}
#endif
#endif
#endif
