/*!
  \file oned_stiffened_fluxes.h
  \brief Riemann solvers for 1D stiffened
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup fluxes
*/
#ifndef ONED_STIFFENED_FLUXES_H
#define ONED_STIFFENED_FLUXES_H
#ifdef ONED
#include <scalar_def.h>
#include <math.h>
#include <macros.h>
#include <basic_fluxes.h>

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
arch_device void oned_stiffened_rusanov(scalar rhoL,
					scalar rhoR,
					scalar vxL,
					scalar vxR,
					scalar EtL,
					scalar EtR,
					scalar alphaL,
					scalar alphaR,
					scalar betaL,
					scalar betaR,
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)                     scalar YL(x), scalar YR(x),
#include "loop.h"
					 scalar nx,
					 scalar* F, scalar* ncterm){

  scalar gammaL = 1.0+1.0/alphaL;
  scalar gammaR = 1.0+1.0/alphaR;
  scalar pinfL = (1-1.0/gammaL)*betaL;
  scalar pinfR = (1-1.0/gammaR)*betaR;
  scalar pL = (gammaL-1)*(EtL - betaL - 0.5*rhoL*vxL*vxL);
  scalar pR = (gammaR-1)*(EtR - betaR - 0.5*rhoR*vxR*vxR);
  scalar aL = sqrt((gammaL*(pL+pinfL))/rhoL);
  scalar aR = sqrt((gammaR*(pR+pinfR))/rhoR);
  int fcnt = 0;   // field counter
  
  // Find the maximum eigenvalue
  scalar maxvap = MAX(fabs(vxL)+aL, fabs(vxR)+aR);
  
  //first: fx = rho*u; 
  F[fcnt] = 0.5*((flux_ab(rhoL,vxL) + flux_ab(rhoR,vxR))*nx
		 -maxvap*(rhoR-rhoL)); fcnt++;
  
  //second: fx = rho*u*u+Bx*Bx+Pbar; 
  F[fcnt] = 0.5*((flux_ab2pc(rhoL,vxL,pL)  + flux_ab2pc(rhoR,vxR,pR))*nx
		 -maxvap*(rhoR*vxR-rhoL*vxL)); fcnt++;
  
  //third: fx = EtplusP*u; 
  F[fcnt] = 0.5*((flux_ab(EtL+pL,vxL) + flux_ab(EtR+pR,vxR))*nx
		 -maxvap*(EtR-EtL)); fcnt++;
  
  //fourth 
  F[fcnt] = 0.5*(-maxvap*(alphaR-alphaL));
  ncterm[fcnt] = -0.5*0.5*(vxL+vxR)*(alphaR-alphaL)*nx; fcnt++;

 
  // pinf
  F[fcnt] = 0.5*(-maxvap*(betaR-betaL));
  ncterm[fcnt] = -0.5*0.5*(vxL+vxR)*(betaR-betaL)*nx; fcnt++;
  
    //mass fractions
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) F[fcnt+x] = 0.5*(flux_ab(YL(x),vxL) + flux_ab(YR(x),vxR))*nx \
    -maxvap*(YR(x)-YL(x)); fcnt++;
#include "loop.h"
  
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
arch_device void oned_stiffened_hll(scalar rhoL,
				    scalar rhoR,
				    scalar vxL,
				    scalar vxR,
				    scalar EtL,
				    scalar EtR,
				    scalar alphaL,
				    scalar alphaR,
				    scalar betaL,
				    scalar betaR,
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)                    scalar YL(x), scalar YR(x),
#include "loop.h"
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
arch_device void oned_stiffened_roe(scalar rhoL,
				    scalar rhoR,
				    scalar vxL,
				    scalar vxR,
				    scalar EtL,
				    scalar EtR,
				    scalar alphaL,
				    scalar alphaR,
				    scalar betaL,
				    scalar betaR,
				    
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)                    scalar YL(x), scalar YR(x),
#include "loop.h"
				    scalar nx,
				    scalar* F, scalar* ncterm){

  scalar gammaL = 1.0+1.0/alphaL;
  scalar gammaR = 1.0+1.0/alphaR;
  scalar pinfL = (1-1.0/gammaL)*betaL;
  scalar pinfR = (1-1.0/gammaR)*betaR;
  scalar pL = (gammaL-1)*(EtL - betaL - 0.5*rhoL*vxL*vxL);
  scalar pR = (gammaR-1)*(EtR - betaR - 0.5*rhoR*vxR*vxR);
  scalar aL = sqrt((gammaL*(pL+pinfL))/rhoL);
  scalar aR = sqrt((gammaR*(pR+pinfR))/rhoR);
  int fcnt = 0; // field counter
  
  // Compute Roe averages
  scalar RT    = sqrt(rhoR/rhoL);
  scalar rho   = RT*rhoL;
  scalar v     = (vxL+RT*vxR)/(1+RT);
  scalar HL    = (EtL + pL)/rhoL;
  scalar HR    = (EtR + pR)/rhoR;
  scalar H     = (HL+RT*HR)/(1+RT);
  scalar alpha = (alphaL+RT*alphaR)/(1+RT);
  scalar gamma = 1+1.0/alpha;
  scalar a     = sqrt((gamma-1)*(H-0.5*v*v));
  scalar iL    = pL*alphaL;
  scalar iR    = pR*alphaR;
  scalar i     = (sqrt(rhoL)*iL+sqrt(rhoR)*iR)/(sqrt(rhoL)+sqrt(rhoR));
  scalar Dp    = (gamma-1)*(gamma-1)*(alpha*(iR-iL) - i*(alphaR-alphaL));
  scalar p     = (gamma-1)*i;

  // Roe waves strengths
  scalar dV0 = (Dp - rho*a*(vxR-vxL))/(2*a*a);
  scalar dV1 = (rhoR-rhoL) - Dp/(a*a);
  scalar dV2 = (Dp + rho*a*(vxR-vxL))/(2*a*a);
  scalar dV3 = alphaR-alphaL;
  scalar dV4 = betaR-betaL;
  
  // Absolute value of Roe eigenvalues
  scalar ws0 = fabs(v-a);
  scalar ws1 = fabs(v);
  scalar ws2 = fabs(v+a);
  scalar ws3 = fabs(v);
  scalar ws4 = fabs(v);
  
  // Entropy fix from http://www.cfdbooks.com/cfdcodes/oned_euler_fluxes_v4.f90
  // Modified wave speeds for nonlinear fields (to remove expansion shocks).
  // There are various ways to implement an entropy fix. This is just one
  // example.
  /* scalar Da = MAX(0.0, 4*((vxR-aR)-(vxL-aL))); */
  /* if(ws[0] < 0.5*Da) ws[0] = ws[0]*ws[0]/Da + 0.25*Da; */
  /* Da = MAX(0.0, 4*((vxR+aR)-(vxL+aL))); */
  /* if(ws[2] < 0.5*Da) ws[2] = ws[2]*ws[2]/Da + 0.25*Da; */

  // Roe Right eigenvectors
  scalar R00 = 1;
  scalar R01 = v-a;
  scalar R02 = H-v*a;
  scalar R03 = 0;
  scalar R04 = 0;

  scalar R10 = 1;
  scalar R11 = v;
  scalar R12 = 0.5*v*v;
  scalar R13 = 0;
  scalar R14 = 0;

  scalar R20 = 1;
  scalar R21 = v+a;
  scalar R22 = H+v*a;
  scalar R23 = 0;
  scalar R24 = 0;

  scalar R30 = 0;
  scalar R31 = 0;
  scalar R32 = p;
  scalar R33 = 1;
  scalar R34 = 0;  

  scalar R40 = 0;
  scalar R41 = 0;
  scalar R42 = 1;
  scalar R43 = 0;
  scalar R44 = 1;  

  // first: fx = rho*u
  F[fcnt] = 0.5*(flux_ab(rhoL,vxL) + flux_ab(rhoR,vxR))*nx
    -0.5*(ws0*dV0*R00+
	  ws1*dV1*R10+
	  ws2*dV2*R20+
	  ws3*dV3*R30+
	  ws4*dV4*R40); fcnt++;

  //second: fx = rho*u*u+Bx*Bx+Pbar; 
  F[fcnt] = 0.5*(flux_ab2pc(rhoL,vxL,pL)  + flux_ab2pc(rhoR,vxR,pR))*nx
    -0.5*(ws0*dV0*R01+
	  ws1*dV1*R11+
	  ws2*dV2*R21+
	  ws3*dV3*R31+
	  ws4*dV4*R41); fcnt++;

  //third: fx = EtplusP*u; 
  F[fcnt] = 0.5*(flux_ab(EtL+pL,vxL) + flux_ab(EtR+pR,vxR))*nx
    -0.5*(ws0*dV0*R02+
	  ws1*dV1*R12+
	  ws2*dV2*R22+
	  ws3*dV3*R32+
	  ws4*dV4*R42); fcnt++;

  //fourth 
  ncterm[fcnt] = -0.5*v*(alphaR-alphaL)*nx;
  F[fcnt] += -0.5*(ws0*dV0*R03+
		   ws1*dV1*R13+
		   ws2*dV2*R23+
		   ws3*dV3*R33+
		   ws4*dV4*R43); fcnt++;

  //fifth: pinf term 
  ncterm[fcnt] = -0.5*v*(betaR-betaL)*nx;
  F[fcnt] += -0.5*(ws0*dV0*R04+
		   ws1*dV1*R14+
		   ws2*dV2*R24+
		   ws3*dV3*R34+
		   ws4*dV4*R44); fcnt++;

  //mass fractions
  scalar Y,dVY; // will hold the Roe average of the mass fraction and the wave strength
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) Y = (YL(x)/rhoL+RT*YR(x)/rhoR)/(1+RT);			\
  dVY = YR(x) - YL(x) - (dV0+dV1+dV2)*Y;				\
  F[fcnt+x] = 0.5*(flux_ab(YL(x),vxL) + flux_ab(YR(x),vxR))*nx		\
    -0.5*(ws0*dV0*Y+							\
	  ws1*dV1*Y+							\
	  ws1*dVY*1+							\
	  ws2*dV2*Y); fcnt++;
#include "loop.h"
  
/*     //mass fractions N-C form */
/* #include "loopstart.h" */
/* #define LOOP_END N_Y */
/* #define MACRO(x) F[4+x] = -0.5*(ws0*dV0*R00+				\ */
/* 				ws1*dV1*R10+				\ */
/* 				ws2*dV2*R20+				\ */
/* 				ws3*dV3*R30);				\ */
/*   ncterm[4+x] = -0.5*v*(YR(x) - YL(x))*nx; */
/* #include "loop.h" */

} // end Roe function
#endif
#endif
#endif
