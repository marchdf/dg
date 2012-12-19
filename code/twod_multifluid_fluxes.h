#ifndef TWOD_MULTIFLUID_FLUXES_H
#define TWOD_MULTIFLUID_FLUXES_H
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
arch_device scalar twod_multifluid_rusanov(scalar* uL,scalar* uR, scalar* n,scalar* F, scalar* ncterm){

  scalar nx = n[0];
  scalar ny = n[1];
  scalar rhoL  = uL[0];
  scalar rhoR  = uR[0];
  scalar vxL   = uL[1]/uL[0];
  scalar vxR   = uR[1]/uR[0];
  scalar vyL   = uL[2]/uL[0];
  scalar vyR   = uR[2]/uR[0];
  scalar vnL = vxL*nx+vyL*ny;
  scalar vnR = vxR*nx+vyR*ny;
  scalar EtL   = uL[3];
  scalar EtR   = uR[3];
#ifdef GAMCONS
  scalar alphaL = uL[4]/uL[0];
  scalar alphaR = uR[4]/uR[0];
#elif  GAMNCON
  scalar alphaL = uL[4];
  scalar alphaR = uR[4];
#endif
  scalar gammaL = 1.0+1.0/alphaL;
  scalar gammaR = 1.0+1.0/alphaR;
  scalar pL = (gammaL-1)*(EtL - 0.5*rhoL*(vxL*vxL+vyL*vyL));
  scalar pR = (gammaR-1)*(EtR - 0.5*rhoR*(vxR*vxR+vyR*vyR));
  scalar aL = sqrt((gammaL*pL)/rhoL);
  scalar aR = sqrt((gammaR*pR)/rhoR);
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
  //fifth
#ifdef GAMCONS
  F[4] = 0.5*(flux_abc(rhoL,vnL,alphaL) + flux_abc(rhoR,vnR,alphaR)
	      -maxvap*(rhoR*alphaR-rhoL*alphaL));
#elif  GAMNCON
  F[4] = 0.5*(-maxvap*(alphaR-alphaL));
  ncterm[4] = -0.5*0.5*(vnL+vnR)*(alphaR-alphaL);
#endif

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
arch_device scalar twod_multifluid_hll(scalar* uL,scalar* uR, scalar* n,scalar* F, scalar* ncterm){

  printf("Not implemented... The code is going to crash.\n");

}

//*****************************************************************************
//* -- Roe's Flux Function ---
//*
//* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
//* Schemes, Journal of Computational Physics, 43, pp. 357-372.
//*
//*****************************************************************************
arch_device scalar twod_multifluid_roe(scalar* uL,scalar* uR, scalar* n,scalar* F, scalar* ncterm){

  scalar nx = n[0];
  scalar ny = n[1];
  scalar tx = -ny;
  scalar ty =  nx;
  scalar rhoL  = uL[0];
  scalar rhoR  = uR[0];
  scalar vxL   = uL[1]/uL[0];
  scalar vxR   = uR[1]/uR[0];
  scalar vyL   = uL[2]/uL[0];
  scalar vyR   = uR[2]/uR[0];
  scalar vnL   = vxL*nx+vyL*ny;
  scalar vnR   = vxR*nx+vyR*ny;
  scalar vtL   = vxL*tx+vyL*ty;
  scalar vtR   = vxR*tx+vyR*ty;
  scalar EtL   = uL[3];
  scalar EtR   = uR[3];
#ifdef GAMCONS
  scalar alphaL = uL[4]/uL[0];
  scalar alphaR = uR[4]/uR[0];
#elif  GAMNCON
  scalar alphaL = uL[4];
  scalar alphaR = uR[4];
#endif
  scalar gammaL = 1.0+1.0/alphaL;
  scalar gammaR = 1.0+1.0/alphaR;
  scalar pL = (gammaL-1)*(EtL - 0.5*rhoL*(vxL*vxL+vyL*vyL));
  scalar pR = (gammaR-1)*(EtR - 0.5*rhoR*(vxR*vxR+vyR*vyR));
  scalar aL = sqrt((gammaL*pL)/rhoL);
  scalar aR = sqrt((gammaR*pR)/rhoR);
  scalar HL = (EtL + pL)/rhoL;
  scalar HR = (EtR + pR)/rhoR;
  scalar iL = pL*alphaL;
  scalar iR = pR*alphaR;

  // Compute Roe averages
  scalar RT  = sqrt(rhoR/rhoL);
  scalar rho = RT*rhoL;
  scalar vx  = (vxL+RT*vxR)/(1+RT);
  scalar vy  = (vyL+RT*vyR)/(1+RT);
  scalar vn  = vx*nx+vy*ny; 
  scalar vt  = vx*tx+vy*ty; 
  scalar H   = (HL+RT*HR)/(1+RT);
  scalar alpha = (alphaL+RT*alphaR)/(1+RT);
  scalar gamma = 1+1.0/alpha;
  scalar a   = sqrt((gamma-1)*(H-0.5*(vx*vx+vy*vy)));
  scalar i   =  (iL+RT*iR)/(1+RT);
  scalar p   = (gamma-1)*i;
  
  // Roe waves strengths
  scalar drho = rhoR - rhoL;
  scalar dp   = (gamma-1)*(gamma-1)*(alpha*(iR-iL) - i*(alphaR-alphaL));
  scalar dvn =  vnR - vnL;
  scalar dvt =  vtR - vtL;
  int sizevap = 5;
  scalar* dV = new scalar[sizevap];
  dV[0] = (dp - rho*a*dvn)/(2*a*a);
  dV[1] = rho*dvt/a;
  dV[2] = drho - dp/(a*a);
  dV[3] = (dp + rho*a*dvn )/(2*a*a);
  dV[4] = alphaR-alphaL;

  // Absolute value of Roe eigenvalues
  scalar* ws = new scalar[sizevap];
  ws[0] = fabs(vn-a);
  ws[1] = fabs(vn);
  ws[2] = fabs(vn);
  ws[3] = fabs(vn+a);
  ws[4] = fabs(vn);

  /* // Entropy fix */


  // Roe Right eigenvectors
  scalar* R = new scalar[sizevap*sizevap];
  R[0*sizevap+0] = 1;
  R[0*sizevap+1] = vx - a*nx;
  R[0*sizevap+2] = vy - a*ny;
  R[0*sizevap+3] = H - vn*a;
  R[0*sizevap+4] = 0;
    
  R[1*sizevap+0] = 0;
  R[1*sizevap+1] = a*tx;
  R[1*sizevap+2] = a*ty;
  R[1*sizevap+3] = vt*a;
  R[1*sizevap+4] = 0;
    
  R[2*sizevap+0] = 1;
  R[2*sizevap+1] = vx;
  R[2*sizevap+2] = vy;
  R[2*sizevap+3] = 0.5*(vx*vx+vy*vy);
  R[2*sizevap+4] = 0;
  
  R[3*sizevap+0] = 1;
  R[3*sizevap+1] = vx + a*nx;
  R[3*sizevap+2] = vy + a*ny;
  R[3*sizevap+3] = H + vn*a;
  R[3*sizevap+4] = 0;

  R[4*sizevap+0] = 0;
  R[4*sizevap+1] = 0;
  R[4*sizevap+2] = 0;
  R[4*sizevap+3] = p;
  R[4*sizevap+4] = 1;

  //first: fx = rho*u; fy = rho*v
  F[0] = 0.5*(flux_ab(rhoL,vnL) + flux_ab(rhoR,vnR));
  for(int k=0;k<sizevap;k++) F[0] += -0.5*ws[k]*dV[k]*R[k*sizevap+0];

  //second: fx = rho*u*u+p; fy = rho*u*v
  F[1] = 0.5*(flux_apb(rhoL*vnL*vxL,pL*nx)  + flux_apb(rhoR*vnR*vxR,pR*nx));
  for(int k=0;k<sizevap;k++) F[1] += -0.5*ws[k]*dV[k]*R[k*sizevap+1];

  //third: fx = rho*u*v; fy = rho*v*v+p
  F[2] = 0.5*(flux_apb(rhoL*vnL*vyL,pL*ny)  + flux_apb(rhoR*vnR*vyR,pR*ny));
  for(int k=0;k<sizevap;k++) F[2] += -0.5*ws[k]*dV[k]*R[k*sizevap+2];
 
  //fourth: fx = rho*u*H; fy = rho*v*H;
  F[3] = 0.5*(flux_abc(rhoL,vnL,HL)+flux_abc(rhoR,vnR,HR));
  for(int k=0;k<sizevap;k++) F[3] += -0.5*ws[k]*dV[k]*R[k*sizevap+3];
  
  //fourth 
#ifdef GAMCONS
  F[4] = 0.5*(flux_abc(rhoL,vnL,alphaL) + flux_abc(rhoR,vnR,alphaR));
#elif  GAMNCON
  F[4] = 0;
  ncterm[4] = -0.5*vn*(alphaR-alphaL);
#endif
  for(int k=0;k<sizevap;k++) F[4] += -0.5*ws[k]*dV[k]*R[k*sizevap+4];

  // Free some pointers
  delete[] dV;
  delete[] ws;
  delete[] R;
  
} // end Roe function


#endif
