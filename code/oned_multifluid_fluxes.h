/*!
  \file oned_multifluid_fluxes.h
  \brief Riemann solvers for 1D multifluid
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \defgroup fluxes Fluxes and Riemann solvers
  \ingroup fluxes
*/
#ifndef ONED_MULTIFLUID_FLUXES_H
#define ONED_MULTIFLUID_FLUXES_H
#ifdef ONED
#include <scalar_def.h>
#include <math.h>
#include <macros.h>
#include <basic_fluxes.h>
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
arch_device void oned_multifluid_rusanov(scalar rhoL,
					 scalar rhoR,
					 scalar vxL,
					 scalar vxR,
					 scalar EtL,
					 scalar EtR,
					 scalar alphaL,
					 scalar alphaR,
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)                     scalar YL(x), scalar YR(x),
#include "loop.h"
					 scalar nx,
					 //scalar x,
					 scalar* F, scalar* ncterm){

#ifdef GAMCONS
  alphaL = alphaL/rhoL;
  alphaR = alphaR/rhoR;
#endif
  scalar gammaL = 1.0+1.0/alphaL;
  scalar gammaR = 1.0+1.0/alphaR;
  scalar pL = (gammaL-1)*(EtL - 0.5*rhoL*vxL*vxL);
  scalar pR = (gammaR-1)*(EtR - 0.5*rhoR*vxR*vxR);
  // Adjust pressure for gravity (flux method)
  //pL = pL - rhoL*constants::GLOBAL_GX*x; //pL = 1e5;
  //pR = pR - rhoR*constants::GLOBAL_GX*x; //pR = 1e5;
  //printf("xf=%f, pL=%f, pR=%f\n",x,pL,pR);  

  scalar aL = sqrt((gammaL*pL)/rhoL);
  scalar aR = sqrt((gammaR*pR)/rhoR);
 
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

  //fourth 
#ifdef GAMCONS
  F[3] = 0.5*((flux_abc(rhoL,vxL,alphaL) + flux_abc(rhoR,vxR,alphaR))*nx
	       -maxvap*(rhoR*alphaR-rhoL*alphaL));
#elif  GAMNCON
  F[3] = 0.5*(-maxvap*(alphaR-alphaL));
  ncterm[3] = -0.5*0.5*(vxL+vxR)*(alphaR-alphaL)*nx;
#endif

    //mass fractions
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) F[4+x] = 0.5*(flux_ab(YL(x),vxL) + flux_ab(YR(x),vxR))*nx \
    -maxvap*(YR(x)-YL(x));
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
arch_device void oned_multifluid_hll(scalar rhoL,
				     scalar rhoR,
				     scalar vxL,
				     scalar vxR,
				     scalar EtL,
				     scalar EtR,
				     scalar alphaL,
				     scalar alphaR,
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)                     scalar YL(x), scalar YR(x),
#include "loop.h"
				     scalar nx,
				     scalar* F, scalar* ncterm){
#ifdef GAMCONS
  alphaL = alphaL/rhoL;
  alphaR = alphaR/rhoR;
#endif
  scalar gammaL = 1.0+1.0/alphaL;
  scalar gammaR = 1.0+1.0/alphaR;
  scalar pL = (gammaL-1)*(EtL - 0.5*rhoL*vxL*vxL);
  scalar pR = (gammaR-1)*(EtR - 0.5*rhoR*vxR*vxR);
  scalar aL = sqrt((gammaL*pL)/rhoL);
  scalar aR = sqrt((gammaR*pR)/rhoR);
  
  // Wave estimates
  scalar SL = 0;
  scalar SR = 0;
  int estimate = 0;
  switch (estimate){
  case 1:{ // Toro 10.49
    scalar RT    = sqrt(rhoR/rhoL);
    scalar v     = (vxL+RT*vxR)/(1+RT);
    scalar HL    = (EtL + pL)/rhoL;
    scalar HR    = (EtR + pR)/rhoR;
    scalar H     = (HL+RT*HR)/(1+RT);
    scalar alpha = (alphaL+RT*alphaR)/(1+RT);
    scalar gamma = 1+1.0/alpha;
    scalar a     = sqrt((gamma-1)*(H-0.5*v*v));
    SL = v*nx-a;
    SR = v*nx+a;
    break;}
  case 2: { // Toro 10.52
    scalar eta2 = 0.5*sqrt(rhoL*rhoR)/((sqrt(rhoL)+sqrt(rhoR))*(sqrt(rhoL)+sqrt(rhoR)));
    scalar dbar2 = (sqrt(rhoL)*aL*aL + sqrt(rhoR)*aR*aR)/(sqrt(rhoL)+sqrt(rhoR)) + eta2*(vxR-vxL)*(vxR-vxL);
    scalar ubar = 0.5*(vxR+vxL);
    SL = ubar*nx - sqrt(dbar2);
    SR = ubar*nx + sqrt(dbar2);
    break;}
  default: // Toro 10.48
    SL = MIN(vxL*nx-aL,vxR*nx-aR); // add 0 to min (see Saurel SIAM 99)
    SR = MAX(vxL*nx+aL,vxR*nx+aR); // add 0 to max (see Saurel SIAM 99)
  }

  // Non-conservative terms (see paper by Rhebergen)
  scalar vnc    = -0.5*    (vxL*nx+vxR*nx)*(alphaL-alphaR);
  scalar vncabs = -0.5*fabs(vxL*nx+vxR*nx)*(alphaL-alphaR);

  // define the flux
  if (SL > 0){
    F[0] = flux_ab(rhoL,vxL)*nx;
    F[1] = flux_ab2pc(rhoL,vxL,pL)*nx;
    F[2] = flux_ab(EtL+pL,vxL)*nx;
#ifdef GAMCONS
    F[3] = flux_abc(rhoL,vxL,alphaL)*nx;
#elif  GAMNCON
    F[3] = -0.5*vncabs;
#endif
    //mass fractions
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) F[4+x] = flux_ab(YL(x),vxL)*nx;
#include "loop.h"
  }
  else if ((SL < 0)&&(SR > 0)){
    F[0] = fhll(rhoL, SL, flux_ab(rhoL,vxL)*nx,   rhoR, SR, flux_ab(rhoR,vxR)*nx);
    F[1] = fhll(  vxL, SL, flux_ab2pc(rhoL,vxL,pL)*nx,  vxR, SR, flux_ab2pc(rhoR,vxR,pR)*nx);
    F[2] = fhll( EtL, SL, flux_ab(EtL+pL,vxL)*nx,    EtR, SR, flux_ab(EtR+pR,vxR)*nx);
#ifdef GAMCONS
    F[3] = fhll(alphaL, SL, flux_abc(rhoL,vxL,alphaL)*nx, alphaR, SR, flux_abc(rhoR,vxR,alphaR)*nx);
#elif  GAMNCON
    F[3] = fhll(alphaL, SL, 0, alphaR, SR, 0) - 0.5*fabs(SR+SL)/fabs(SR-SL)*vncabs;
#endif
    //mass fractions
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) F[4+x] = fhll(YL(x), SL, flux_ab(YL(x),vxL)*nx,   YR(x), SR, flux_ab(YR(x),vxR)*nx);
#include "loop.h"
  }
  else if (SR < 0){
    F[0] = flux_ab(rhoR,vxR)*nx;
    F[1] = flux_ab2pc(rhoR,vxR,pR)*nx;
    F[2] = flux_ab(EtR+pR,vxR)*nx;
#ifdef GAMCONS
    F[3] = flux_abc(rhoR,vxR,alphaR)*nx;
#elif  GAMNCON
    F[3] = 0.5*vncabs;
#endif
    //mass fractions
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) F[4+x] = flux_ab(YR(x),vxR)*nx;
#include "loop.h"
  }

#ifdef GAMNCON
  ncterm[3] = -0.5*vnc;
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
arch_device void oned_multifluid_roe(scalar rhoL,
				     scalar rhoR,
				     scalar vxL,
				     scalar vxR,
				     scalar EtL,
				     scalar EtR,
				     scalar alphaL,
				     scalar alphaR,
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)                     scalar YL(x), scalar YR(x),
#include "loop.h"
				     scalar nx,
				     scalar* F, scalar* ncterm){

#ifdef GAMCONS
  alphaL = alphaL/rhoL;
  alphaR = alphaR/rhoR;
#endif
  scalar gammaL = 1.0+1.0/alphaL;
  scalar gammaR = 1.0+1.0/alphaR;
  scalar pL = (gammaL-1)*(EtL - 0.5*rhoL*vxL*vxL);
  scalar pR = (gammaR-1)*(EtR - 0.5*rhoR*vxR*vxR);
  // Adjust pressure for gravity (flux method)
  //pL = pL - rhoL*constants::GLOBAL_GX*x; //pL = 1e5;
  //pR = pR - rhoR*constants::GLOBAL_GX*x; //pR = 1e5;
  //printf("xf=%f, pL=%f, pR=%f\n",x,pL,pR);  
  
  scalar aL = sqrt((gammaL*pL)/rhoL);
  scalar aR = sqrt((gammaR*pR)/rhoR);


  
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
  
  // Absolute value of Roe eigenvalues
  scalar ws0 = fabs(v-a);
  scalar ws1 = fabs(v);
  scalar ws2 = fabs(v+a);
  scalar ws3 = fabs(v);
  
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

  scalar R10 = 1;
  scalar R11 = v;
  scalar R12 = 0.5*v*v;
  scalar R13 = 0;

  scalar R20 = 1;
  scalar R21 = v+a;
  scalar R22 = H+v*a;
  scalar R23 = 0;

  scalar R30 = 0;
  scalar R31 = 0;
  scalar R32 = p;
  scalar R33 = 1;  

  // first: fx = rho*u
  F[0] = 0.5*(flux_ab(rhoL,vxL) + flux_ab(rhoR,vxR))*nx
    -0.5*(ws0*dV0*R00+
	  ws1*dV1*R10+
	  ws2*dV2*R20+
	  ws3*dV3*R30);

  //second: fx = rho*u*u+Bx*Bx+Pbar; 
  F[1] = 0.5*(flux_ab2pc(rhoL,vxL,pL)  + flux_ab2pc(rhoR,vxR,pR))*nx
    -0.5*(ws0*dV0*R01+
	  ws1*dV1*R11+
	  ws2*dV2*R21+
	  ws3*dV3*R31);

  //third: fx = EtplusP*u; 
  F[2] = 0.5*(flux_ab(EtL+pL,vxL) + flux_ab(EtR+pR,vxR))*nx
    -0.5*(ws0*dV0*R02+
	  ws1*dV1*R12+
	  ws2*dV2*R22+
	  ws3*dV3*R32);

  //fourth 
#ifdef GAMCONS
  F[3] = 0.5*(flux_abc(rhoL,vxL,alphaL) + flux_abc(rhoR,vxR,alphaR))*nx;
#elif  GAMNCON
  F[3] = 0;
  ncterm[3] = -0.5*v*(alphaR-alphaL)*nx;
#endif
  F[3] += -0.5*(ws0*dV0*R03+
		ws1*dV1*R13+
		ws2*dV2*R23+
		ws3*dV3*R33);

  //mass fractions
  scalar Y,dVY; // will hold the Roe average of the mass fraction and the wave strength
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) Y = (YL(x)/rhoL+RT*YR(x)/rhoR)/(1+RT);			\
  dVY = YR(x) - YL(x) - (dV0+dV2)*Y;					\
  F[4+x] = 0.5*(flux_ab(YL(x),vxL) + flux_ab(YR(x),vxR))*nx		\
    -0.5*(ws0*dV0*Y+							\
	  ws1*dVY*1+							\
	  ws2*dV2*Y);
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
