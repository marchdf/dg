#ifndef TWOD_PASSIVE_FLUXES_H
#define TWOD_PASSIVE_FLUXES_H
#include <scalar_def.h>
#include <math.h>
#include <macros.h>
#include <constants.h>
#include <basic_fluxes.h>


//*****************************************************************************
//* --- Rusanov's Flux Function ---
//*
//* V. V. Rusanov, Calculation of Interaction of Non-Steady Shock Waves with
//* Obstacles, J. Comput. Math. Phys. USSR, 1, pp. 267-279, 1961.
//*
//*****************************************************************************
arch_device scalar twod_passive_rusanov(scalar* uL,scalar* uR, scalar* n,scalar* F, scalar* ncterm){

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
  scalar phicL = uL[4]/uL[0];
  scalar phicR = uR[4]/uR[0];
  scalar phincL= uL[5];
  scalar phincR= uR[5];
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

  //fifth: fx = rho*u*phi
  F[4] = 0.5*(flux_abc(rhoL,vnL,phicL) + flux_abc(rhoR,vnR,phicR)
	      -maxvap*(rhoR*phicR-rhoL*phicL));
	
  //sixth:
  F[5] = 0.5*(-maxvap*(phincR-phincL));
  ncterm[5] = 0.5*0.5*(vnL+vnR)*(phincL-phincR);

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
arch_device scalar twod_passive_hll(scalar* uL,scalar* uR, scalar* n,scalar* F, scalar* ncterm){

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
  scalar phicL = uL[4]/uL[0];
  scalar phicR = uR[4]/uR[0];
  scalar phincL= uL[5];
  scalar phincR= uR[5];
  scalar gamma = constants::GLOBAL_GAMMA;
  scalar pL = (gamma-1)*(EtL - 0.5*rhoL*(vxL*vxL+vyL*vyL));
  scalar pR = (gamma-1)*(EtR - 0.5*rhoR*(vxR*vxR+vyR*vyR));
  scalar aL = sqrt((gamma*pL)/rhoL);
  scalar aR = sqrt((gamma*pR)/rhoR);

  // Wave estimates
  scalar SL = 0;
  scalar SR = 0;
  int estimate = 0;
  switch (estimate){
  case 1:{ // Toro 10.49
    scalar RT   = sqrt(rhoR/rhoL);
    scalar uRoe = (vnL+RT*vnR)/(1+RT);
    scalar HL   = (EtL + pL)/rhoL;
    scalar HR   = (EtR + pR)/rhoR;
    scalar HRoe = (HL+RT*HR)/(1+RT);
    scalar aRoe = sqrt((gamma-1)*(HRoe-0.5*uRoe*uRoe));
    SL = uRoe-aRoe;
    SR = uRoe+aRoe;
    break;}
  case 2: { // Toro 10.52
    scalar eta2 = 0.5*sqrt(rhoL*rhoR)/((sqrt(rhoL)+sqrt(rhoR))*(sqrt(rhoL)+sqrt(rhoR)));
    scalar dbar2 = (sqrt(rhoL)*aL*aL + sqrt(rhoR)*aR*aR)/(sqrt(rhoL)+sqrt(rhoR)) + eta2*(vxR-vxL)*(vxR-vxL);
    scalar ubar = 0.5*(vnR+vnL);
    SL = ubar - sqrt(dbar2);
    SR = ubar + sqrt(dbar2);
    break;}
  default: // Toro 10.48
    SL = MIN(vnL-aL,vnR-aR);
    SR = MAX(vnL+aL,vnR+aR);
  }

  // Non-conservative terms (see paper by Rhebergen)
  scalar vnc    = -0.5*    (vnL+vnR)*(phincL-phincR);
  scalar vncabs = -0.5*fabs(vnL+vnR)*(phincL-phincR);

  // define the flux
  if (SL > 0){
    F[0] = flux_ab(rhoL,vxL)       *nx + flux_ab(rhoL,vyL)       *ny;
    F[1] = flux_ab2pc(rhoL,vxL,pL) *nx + flux_abc(rhoL,vxL,vyL)  *ny;
    F[2] = flux_abc(rhoL,vxL,vyL)  *nx + flux_ab2pc(rhoL,vyL,pL) *ny;
    F[3] = flux_ab(EtL+pL,vxL)     *nx + flux_ab(EtL+pL,vyL)     *ny;
    F[4] = flux_abc(rhoL,vxL,phicL)*nx + flux_abc(rhoL,vyL,phicL)*ny;
    F[5] = -0.5*vncabs;
  }
  else if ((SL < 0)&&(SR > 0)){
    F[0] = fhll(  rhoL, SL, flux_ab(rhoL,vxL)       *nx + flux_ab(rhoL,vyL)       *ny,   rhoR, SR, flux_ab(rhoR,vxR)       *nx + flux_ab(rhoR,vyR)       *ny);
    F[1] = fhll(   vxL, SL, flux_ab2pc(rhoL,vxL,pL) *nx + flux_abc(rhoL,vxL,vyL)  *ny,    vxR, SR, flux_ab2pc(rhoR,vxR,pR) *nx + flux_abc(rhoR,vxR,vyR)  *ny);
    F[2] = fhll(   vyL, SL, flux_abc(rhoL,vxL,vyL)  *nx + flux_ab2pc(rhoL,vyL,pL) *ny,    vyR, SR, flux_ab2pc(rhoR,vxR,pR) *nx + flux_ab2pc(rhoR,vyR,pR) *ny);
    F[3] = fhll(   EtL, SL, flux_ab(EtL+pL,vxL)     *nx + flux_ab(EtL+pL,vyL)     *ny,    EtR, SR, flux_ab(EtR+pR,vxR)     *nx + flux_ab(EtR+pR,vyR)     *ny);
    F[4] = fhll( phicL, SL, flux_abc(rhoL,vxL,phicL)*nx + flux_abc(rhoL,vyL,phicL)*ny,  phicR, SR, flux_abc(rhoR,vxR,phicR)*nx + flux_abc(rhoR,vyR,phicR)*ny);
    F[5] = fhll(phincL, SL, 0                                                        , phincR, SR, 0                                                        )
      - 0.5*fabs(SR+SL)/fabs(SR-SL)*vncabs;
  }
  else if (SR < 0){
    F[0] = flux_ab(rhoR,vxR)       *nx  + flux_ab(rhoR,vyR)       *ny;
    F[1] = flux_ab2pc(rhoR,vxR,pR) *nx  + flux_abc(rhoR,vxR,vyR)  *ny;
    F[2] = flux_abc(rhoR,vxR,vyR)  *nx  + flux_ab2pc(rhoR,vyR,pR) *ny;
    F[3] = flux_ab(EtR+pR,vxR)     *nx  + flux_ab(EtR+pR,vyR)     *ny;
    F[4] = flux_abc(rhoR,vxR,phicR)*nx  + flux_abc(rhoR,vyR,phicR)*ny;
    F[5] = 0.5*vncabs;
  }
  ncterm[5] = -0.5*vnc;

} // end HLL function


//*****************************************************************************
//* -- Roe's Flux Function ---
//*
//* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
//* Schemes, Journal of Computational Physics, 43, pp. 357-372.
//*
//*****************************************************************************
arch_device scalar twod_passive_roe(scalar* uL,scalar* uR, scalar* n,scalar* F, scalar* ncterm){

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
  scalar phicL = uL[4]/uL[0];
  scalar phicR = uR[4]/uR[0];
  scalar phincL= uL[5];
  scalar phincR= uR[5];
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
  int sizevap = 4;
  scalar* dV = new scalar[sizevap];
  dV[0] = (dp - rho*a*dvn )/(2*a*a);
  dV[1] = rho*dvt/a;
  dV[2] = drho - dp/(a*a);
  dV[3] = (dp + rho*a*dvn )/(2*a*a);

  // Absolute value of Roe eigenvalues
  scalar* ws = new scalar[sizevap];
  ws[0] = fabs(vn-a);
  ws[1] = fabs(vn);
  ws[2] = fabs(vn);
  ws[3] = fabs(vn+a);

  // Harten's Entropy Fix JCP(1983), 49, pp357-393:
  // only for the nonlinear fields.
  // http://www.cfdbooks.com/cfdcodes/twod_euler_fluxes_v2.f90
  /* scalar dws0 = 0.2; */
  /* if(ws[0] < dws0) ws[0] = 0.5 * (ws[0]*ws[0]/dws0+dws0); */
  /* scalar dws3 = 0.2; */
  /* if(ws[3] < dws3) ws[3] = 0.5 * (ws[3]*ws[3]/dws3+dws3); */

  // Roe Right eigenvectors
  scalar* R = new scalar[sizevap*sizevap];
  R[0*sizevap+0] = 1;
  R[0*sizevap+1] = vx - a*nx;
  R[0*sizevap+2] = vy - a*ny;
  R[0*sizevap+3] = H - vn*a;

  R[1*sizevap+0] = 0;
  R[1*sizevap+1] = a*tx;
  R[1*sizevap+2] = a*ty;
  R[1*sizevap+3] = vt*a;

  R[2*sizevap+0] = 1;
  R[2*sizevap+1] = vx;
  R[2*sizevap+2] = vy;
  R[2*sizevap+3] = 0.5*(vx*vx+vy*vy);

  R[3*sizevap+0] = 1;
  R[3*sizevap+1] = vx + a*nx;
  R[3*sizevap+2] = vy + a*ny;
  R[3*sizevap+3] = H + vn*a;

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
  
  //fifth: fx = rho*u*phi
  F[4] = 0.5*(flux_abc(rhoL,vnL,phicL) + flux_abc(rhoR,vnR,phicR));
  for(int k=0;k<sizevap;k++) F[4] += -0.5*ws[k]*dV[k]*R[k*sizevap+0];
  
  //sixth:
  F[5] = 0.0;
  for(int k=0;k<sizevap;k++) F[5] += -0.5*ws[k]*dV[k]*R[k*sizevap+0];
  ncterm[5] = 0.5*vn*(phincL-phincR);

  // Free some pointers
  delete[] dV;
  delete[] ws;
  delete[] R;
  
} // end Roe function

#endif 
