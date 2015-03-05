/*!
  \file twod_singlefluid_fluxes.h
  \brief Riemann solvers for 2D single fluid Euler
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup fluxes
*/
#ifndef TWOD_SINGLEFLUID_FLUXES_H
#define TWOD_SINGLEFLUID_FLUXES_H
#ifdef TWOD
#include <scalar_def.h>
#include <math.h>
#include <macros.h>
#include <constants.h>
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
arch_device void twod_singlefluid_rusanov(scalar rhoL,
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
arch_device void twod_singlefluid_hll(scalar rhoL,
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

#ifdef USE_CPU
  printf("Not implemented... The code is going to crash.\n");
#endif
  
}

//*****************************************************************************
//* -- Roe's Flux Function ---
//*
//* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
//* Schemes, Journal of Computational Physics, 43, pp. 357-372.
//*
//*****************************************************************************
#elif ROE
arch_device void twod_singlefluid_roe(scalar rhoL,
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
  
} // end Roe function
#endif
#endif
#endif
