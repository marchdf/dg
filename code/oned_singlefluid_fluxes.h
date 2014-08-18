/*!
  \file oned_singlefluid_fluxes.h
  \brief Riemann solvers for 1D Euler equations (singe fluids)
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \defgroup fluxes Fluxes and Riemann solvers
  \ingroup fluxes
*/
#ifndef ONED_SINGLEFLUID_FLUXES_H
#define ONED_SINGLEFLUID_FLUXES_H
#ifdef ONED
#include <scalar_def.h>
#include <math.h>
#include <macros.h>
#include <basic_fluxes.h>
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
#endif
#endif
#endif
