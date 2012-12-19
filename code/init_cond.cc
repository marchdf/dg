#include <init_cond.h>

void init_dg_shallow(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  if (N_F!=3) printf("You are setting up the wrong problem. N_F =%i != 3.\n",N_F);
  
  // Initial conditions
  for(int i = 0; i < N_s; i++){
    for(int e = 0; e < N_E; e++){
      scalar x = XYZNodes(i,e*D+0);
      scalar y = XYZNodes(i,e*D+1);
      scalar sxc= 0.5*0.5;
      scalar syc= 0.5*0.5;
      //U(i,e*N_F+0) = 0.5*(scalar)exp(-(0.5*(x+0.4)*(x+0.4)/sxc) - 0.5*(y+0.4)*(y+0.4)/syc);
      U(i,e*N_F+0) = (scalar)exp(-(0.5*x*x/sxc) - 0.5*y*y/syc);
    }
  }
}


void init_dg_simplew_multifluid(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions (see CFD2 project 1, Sreenivas)

  scalar a = 0;
  scalar u = 0;
  scalar rho = 0;
  scalar gamma = 1.4;
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if (x<=-1.5){
	u = -2.0/gamma;
	a = (1-gamma)/gamma + 1;
	rho = gamma*pow(a,2.0/(gamma-1));
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
      }
      else if ((x>-1.5)&&(x<-0.5)){
	u = -1.0/gamma*(1-tanh((x+1)/(0.25-(x+1)*(x+1))));
	a = u*(gamma-1)/2.0 + 1.0;
	rho = gamma*pow(a,2.0/(gamma-1));;
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
      }
      else if ((x>=-0.5)&&(x<=0.5)){
	u = 0;
	a = 0;
	rho = gamma;
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)1.0/(gamma-1);
      }
      else if ((x>0.5)&&(x<1.5)){
	u = 1.0/gamma*(1+tanh((x-1)/(0.25-(x-1)*(x-1))));;
	a = 1 - (gamma-1)/2.0*u;
	rho = gamma*pow(a,2.0/(gamma-1));
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
      }
      else if (x>=1.5){
	u = 2.0/gamma;
	a = 1 - (gamma-1)/gamma;
	rho = gamma*pow(a,2.0/(gamma-1));
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
      }
#ifdef GAMCONS
      U(i,e*N_F+3) = (scalar)rho/(gamma-1);
#elif GAMNCON
      U(i,e*N_F+3) = (scalar)1.0/(gamma-1);
#endif
    }
  }
}

void buildLRstates_multifluid(scalar rhoL, scalar uL, scalar EtL, scalar gammaL, scalar rhoR, scalar uR, scalar EtR, scalar gammaR, const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  scalar* L = new scalar[N_F];
  scalar* R = new scalar[N_F];

#ifdef ONED
  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);
  L[0] = rhoL;            R[0] = rhoR; 
  L[1] = rhoL*uL;         R[1] = rhoR*uR; 
  L[2] = EtL;             R[2] = EtR;   
#ifdef GAMCONS
  L[3] = rhoL/(gammaL-1); R[3] = rhoR/(gammaR-1);
#elif GAMNCON
  L[3] = 1.0/(gammaL-1);  R[3] = 1.0/(gammaR-1);
#endif
#elif TWOD
  if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);
  scalar vL = 0;          scalar vR = 0;
  EtL  = EtL+ 0.5*rhoL*vL*vL;
  EtR  = EtR+ 0.5*rhoR*vR*vR;
  L[0] = rhoL;            R[0] = rhoR; 
  L[1] = rhoL*uL;         R[1] = rhoR*uR; 
  L[2] = rhoL*vL;         R[2] = rhoL*vL; 
  L[3] = EtL;             R[3] = EtR; 
#ifdef GAMCONS
  L[4] = rhoL/(gammaL-1); R[4] = rhoR/(gammaR-1);
#elif GAMNCON
  L[4] = 1.0/(gammaL-1);  R[4] = 1.0/(gammaR-1);
#endif
  
#endif
  
  for(int e = 0; e < N_E; e++){
    scalar x = XYZNodes(0,e*D+0);
    for(int i = 0; i < N_s; i++){
      if (x<1E-8){
	for(int k = 0; k < N_F; k++) U(i,e*N_F+k) = L[k];
      }
      else if (x>=1E-8){
	for(int k = 0; k < N_F; k++) U(i,e*N_F+k) = R[k];
      }
    }
  }

  delete[] L;
  delete[] R;

}

void init_dg_sodtube_multifluid(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  // Initial conditions
  // U = (  rho, rho ux, rho uy, rho uz,   Bx, By, Bz,    E,   ee)
  //   = (    1,      0,      0,      0,   0,  0,  0, 1.78,  0.5)  for (x<0)
  //   = (0.125,      0,      0,      0,   0,  0,  0, 0.88, 0.05)  for (x>=0)
  // gamma = 1.4

  // Left state
  scalar rhoL = 1;
  scalar uL   = 0;
  scalar pL   = 1.0;
  scalar gammaL= 1.4;
  scalar EtL  = 1.0/(gammaL-1.0)*pL + 0.5*rhoL*uL*uL;
    
  // Right state
  scalar rhoR = 0.125;
  scalar uR   = 0;
  scalar pR   = 0.1;
  scalar gammaR= 1.6;
  scalar EtR  = 1.0/(gammaR-1.0)*pR + 0.5*rhoR*uR*uR;
  
  buildLRstates_multifluid(rhoL, uL, EtL, gammaL, rhoR, uR, EtR, gammaR, N_s, N_E, N_F, D, XYZNodes, U);
}

void init_dg_sodmono_multifluid(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  // Initial conditions
  // U = (  rho, rho ux, rho uy, rho uz,   Bx, By, Bz,    E,   ee)
  //   = (    1,      0,      0,      0,   0,  0,  0, 1.78,  0.5)  for (x<0)
  //   = (0.125,      0,      0,      0,   0,  0,  0, 0.88, 0.05)  for (x>=0)
  // gamma = 1.4

  // Left state
  scalar rhoL = 1;
  scalar uL   = 0;
  scalar pL   = 1.0;
  scalar gammaL= 1.4;
  scalar EtL  = 1.0/(gammaL-1.0)*pL + 0.5*rhoL*uL*uL;
    
  // Right state
  scalar rhoR = 0.125;
  scalar uR   = 0;
  scalar pR   = 0.1;
  scalar gammaR= 1.4;
  scalar EtR  = 1.0/(gammaR-1.0)*pR + 0.5*rhoR*uR*uR;
  
  buildLRstates_multifluid(rhoL, uL, EtL, gammaL, rhoR, uR, EtR, gammaR, N_s, N_E, N_F, D, XYZNodes, U);
}


void init_dg_contact_multifluid(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  // Left state
  scalar rhoL  = 1.0;
  scalar uL    = 1.0;
  scalar gammaL= 1.4;
  scalar pL    = 1.0;
  scalar EtL   = 1.0/(gammaL-1.0)*pL + 0.5*rhoL*uL*uL;
    
  // Right state
  scalar rhoR   = 1.0;
  scalar uR     = 1.0;
  scalar gammaR = 1.4;
  scalar pR     = 1.0;
  scalar EtR    = 1.0/(gammaR-1.0)*pR + 0.5*rhoR*uR*uR;

  buildLRstates_multifluid(rhoL, uL, EtL, gammaL, rhoR, uR, EtR, gammaR, N_s, N_E, N_F, D, XYZNodes, U);

}

void init_dg_rhotact_multifluid(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){
  
  // Left state
  scalar rhoL  = 1.0;
  scalar uL    = 1.0;
  scalar gammaL= 1.4;
  scalar pL    = 1.0;
  scalar EtL   = 1.0/(gammaL-1.0)*pL + 0.5*rhoL*uL*uL;
  
  // Right state
  scalar rhoR   = 0.125;
  scalar uR     = 1.0;
  scalar gammaR = 1.4;
  scalar pR     = 1.0;
  scalar EtR    = 1.0/(gammaR-1.0)*pR + 0.5*rhoR*uR*uR;

  buildLRstates_multifluid(rhoL, uL, EtL, gammaL, rhoR, uR, EtR, gammaR, N_s, N_E, N_F, D, XYZNodes, U);
}


void init_dg_matfrnt_multifluid(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  // Left state
  scalar rhoL  = 1.0;
  scalar uL    = 1.0;
  scalar gammaL= 1.4;
  scalar pL    = 1.0;
  scalar EtL   = 1.0/(gammaL-1.0)*pL + 0.5*rhoL*uL*uL;
    
  // Right state
  scalar rhoR   = 0.125;
  scalar uR     = 1.0;
  scalar gammaR = 1.6;
  scalar pR     = 1.0;
  scalar EtR    = 1.0/(gammaR-1.0)*pR + 0.5*rhoR*uR*uR;

  buildLRstates_multifluid(rhoL, uL, EtL, gammaL, rhoR, uR, EtR, gammaR, N_s, N_E, N_F, D, XYZNodes, U);
  
}

void init_dg_sinegam_multifluid(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  scalar rho     = 1.0;
  scalar u       = 1.0;
  scalar gamma0  = 1.4;
  scalar sinegam = 0.0;  // sine perturbation on gamma
  scalar sinerho = 0.0;  // sine perturbation on rho
  scalar Agam    = 0.20; // amplitude of the perturbation on gamma
  scalar Arho    = 0.20; // amplitude of the perturbation on rho
  scalar p       = 1.0;
  
  scalar* Q = new scalar[N_F];

  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      sinegam = Agam*sin(M_PI*x);
      sinerho = Arho*sin(4*M_PI*x);
#ifdef ONED
      if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);
      Q[0] = rho+sinerho;
      Q[1] = (rho+sinerho)*u;
      Q[2] = 1.0/(gamma0+sinegam-1.0)*p + 0.5*(rho+sinerho)*u*u;
#ifdef GAMCONS
      Q[3] = (scalar)(rho+sinerho)/(gamma0+sinegam-1);
#elif GAMNCON
      Q[3] = (scalar)1.0/(gamma0+sinegam-1);
#endif
#elif TWOD
      if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);
      scalar v = 0;
      Q[0] = rho+sinerho;
      Q[1] = (rho+sinerho)*u;
      Q[2] = (rho+sinerho)*v;
      Q[3] = 1.0/(gamma0+sinegam-1.0)*p + 0.5*(rho+sinerho)*(u*u+v*v);
#ifdef GAMCONS
      Q[4] = (scalar)(rho+sinerho)/(gamma0+sinegam-1);
#elif GAMNCON
      Q[4] = (scalar)1.0/(gamma0+sinegam-1);
#endif
#endif
      for(int k = 0; k < N_F; k++) U(i,e*N_F+k) = Q[k];
    }
  }

  delete[] Q;
}
             
void init_dg_expogam_multifluid(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);

  scalar rho     = 1.0;
  scalar u       = 1.0;
  scalar gamma0  = 1.4;
  scalar expogam = 0.0;  // exp perturbation on gamma
  scalar A       = 0.10; // amplitude of the perturbation
  scalar p       = 1.0;

  scalar* Q = new scalar[N_F];
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      expogam = A*exp(-x*x/(0.2*0.2));
#ifdef ONED
      if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);
      Q[0] = rho;
      Q[1] = rho*u;
      Q[2] = 1.0/(gamma0+expogam-1.0)*p + 0.5*rho*u*u;
#ifdef GAMCONS
      Q[3] = (scalar)rho/(gamma0+expogam-1);
#elif GAMNCON
      Q[3] = (scalar)1.0/(gamma0+expogam-1);
#endif
#elif TWOD
      if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);
      scalar v = 0;
      Q[0] = rho;
      Q[1] = rho*u;
      Q[2] = rho*v;
      Q[3] = 1.0/(gamma0+expogam-1.0)*p + 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
      Q[4] = (scalar)rho/(gamma0+expogam-1);
#elif GAMNCON
      Q[4] = (scalar)1.0/(gamma0+expogam-1);
#endif
#endif
      for(int k = 0; k < N_F; k++) U(i,e*N_F+k) = Q[k];
    }
  }

  delete[] Q;
}

void init_dg_shckint_multifluid(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  // Pre-shock state (material 1)
  scalar rho02   = 0.1;
  scalar u02     =-2.0;
  scalar gamma02 = 1.6667;
  scalar p02     = 1.0;
  scalar Et02    = 1.0/(gamma02-1.0)*p02 + 0.5*rho02*u02*u02;
  scalar c02     = sqrt(gamma02*p02/rho02); // sound speed
  scalar M02     = u02/c02; //Mach number ahead of shock
  
  // Pre-shock state (material 2)
  scalar rho01   = 1.0;
  scalar u01     =-2.0;
  scalar gamma01 = 1.4;
  scalar p01     = 1.0;
  scalar Et01    = 1.0/(gamma01-1.0)*p01 + 0.5*rho01*u01*u01;

  // Post-shock state (material 1) (see p 101 Toro)
  //scalar Ms = 9;   // Shock Mach number
  scalar Ms     = M02+sqrt((gamma02+1)/(2.0*gamma02)*100.0 + (gamma02-1)/(2.0*gamma02));   // Shock Mach number (with ratio)
  printf("Ms=%f\n",Ms);
  scalar s      = Ms*c02;  // shock speed
  scalar rho4   = rho02*(gamma02+1) * (M02-Ms)*(M02-Ms)/((gamma02-1) * (M02-Ms)*(M02-Ms) + 2);
  scalar p4     = p02  *(2*gamma02*(M02-Ms)*(M02-Ms) - (gamma02-1))/(gamma02+1);
  scalar u4     = (1 - rho02/rho4)*s + u02*rho02/rho4;
  scalar gamma4 = gamma02;
  scalar Et4    = 1.0/(gamma4-1.0)*p4 + 0.5*rho4*u4*u4;

  int N_S = 3; // number of states
  fullMatrix<scalar> xS(N_S+1,1); // state boundaries
  xS(0,0) = -1000; xS(1,0) = -0.8; xS(2,0) = -0.2; xS(3,0) = 1000;
  fullMatrix<scalar> S(N_F,N_S); // holds the states
  
#ifdef ONED
  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);
  S(0,0) = rho4;            S(0,1) = rho02;             S(0,2) = rho01;
  S(1,0) = rho4*u4;         S(1,1) = rho02*u02;         S(1,2) = rho01*u01;
  S(2,0) = Et4;             S(2,1) = Et02;              S(2,2) = Et01; 
#ifdef GAMCONS
  S(3,0) = rho4/(gamma4-1); S(3,1) = rho02/(gamma02-1); S(3,2) = rho01/(gamma01-1); 
#elif GAMNCON
  S(3,0) = 1.0/(gamma4-1);  S(3,1) = 1.0/(gamma02-1);   S(3,2) = 1.0/(gamma01-1); 
#endif
#elif TWOD
  if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);
  scalar v4 = 0;            scalar v02 = 0;             scalar v01 = 0;
  Et4  = 1.0/(gamma4-1.0)*p4 + 0.5*rho4*(u4*u4+v4*v4);
  Et02  = 1.0/(gamma02-1.0)*p02 + 0.5*rho02*(u02*u02+v02*v02);
  Et01  = 1.0/(gamma01-1.0)*p01 + 0.5*rho01*(u01*u01+v01*v01);
  S(0,0) = rho4;            S(0,1) = rho02;             S(0,2) = rho01;
  S(1,0) = rho4*u4;         S(1,1) = rho02*u02;         S(1,2) = rho01*u01;
  S(2,0) = rho4*v4;         S(2,1) = rho02*v02;         S(2,2) = rho01*v01;
  S(3,0) = Et4;             S(3,1) = Et02;              S(3,2) = Et01; 
#ifdef GAMCONS
  S(4,0) = rho4/(gamma4-1); S(4,1) = rho02/(gamma02-1); S(4,2) = rho01/(gamma01-1); 
#elif GAMNCON
  S(4,0) = 1.0/(gamma4-1);  S(4,1) = 1.0/(gamma02-1);   S(4,2) = 1.0/(gamma01-1); 
#endif
#endif

  int ind = 0; // index for the state we use
  for(int e = 0; e < N_E; e++){
    scalar x = XYZNodes(0,e*D+0);
    for(int i = 0; i < N_s; i++){
      for(int b = 0; b < N_S; b++) if ((xS(b,0) <= x) && (x < xS(b+1,0)))  ind = b;
      for(int k = 0; k < N_F; k++) U(i,e*N_F+k) = S(k,ind);
    }
  }
}

void init_dg_multint_multifluid(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  // Pre-shock state (material 1)
  scalar rho02   = 0.1;
  scalar u02     =-2.0;
  scalar gamma02 = 1.4;
  scalar p02     = 1.0;
  scalar Et02    = 1.0/(gamma02-1.0)*p02 + 0.5*rho02*u02*u02;
  scalar c02     = sqrt(gamma02*p02/rho02); // sound speed
  scalar M02     = u02/c02; //Mach number ahead of shock
  
  // Pre-shock state (material 2)
  scalar rho01   = 1.0;
  scalar u01     =-2.0;
  scalar gamma01 = 1.1;
  scalar p01     = 1.0;
  scalar Et01    = 1.0/(gamma01-1.0)*p01 + 0.5*rho01*u01*u01;

  // Pre-shock state (material 3)
  scalar rho03   = 1.0;
  scalar u03     =-2.0;
  scalar gamma03 = 1.6667;
  scalar p03     = 1.0;
  scalar Et03    = 1.0/(gamma03-1.0)*p03 + 0.5*rho03*u03*u03;

  // Post-shock state (material 1) (see p 101 Toro)
  //scalar Ms = 9;   // Shock Mach number
  scalar Ms     = M02+sqrt((gamma02+1)/(2.0*gamma02)*100.0 + (gamma02-1)/(2.0*gamma02));   // Shock Mach number (with ratio)
  scalar s      = Ms*c02;  // shock speed
  scalar rho4   = rho02*(gamma02+1) * (M02-Ms)*(M02-Ms)/((gamma02-1) * (M02-Ms)*(M02-Ms) + 2);
  scalar p4     = p02  *(2*gamma02*(M02-Ms)*(M02-Ms) - (gamma02-1))/(gamma02+1);
  scalar u4     = (1 - rho02/rho4)*s + u02*rho02/rho4;
  scalar gamma4 = gamma02;
  scalar Et4    = 1.0/(gamma4-1.0)*p4 + 0.5*rho4*u4*u4;

  int N_S = 4; // number of states
  fullMatrix<scalar> xS(N_S+1,1); // state boundaries
  xS(0,0) = -1000; xS(1,0) = -0.8; xS(2,0) = -0.2; xS(3,0) = 0.3; xS(4,0) = 1000;
  fullMatrix<scalar> S(N_F,N_S); // holds the states
  
#ifdef ONED
  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);
  S(0,0) = rho4;            S(0,1) = rho02;             S(0,2) = rho01;             S(0,3) = rho01;		   
  S(1,0) = rho4*u4;         S(1,1) = rho02*u02;         S(1,2) = rho01*u01;	    S(1,3) = rho01*u01;	   
  S(2,0) = Et4;             S(2,1) = Et02;              S(2,2) = Et01; 		    S(2,3) = Et01; 		   
#ifdef GAMCONS									                               
  S(3,0) = rho4/(gamma4-1); S(3,1) = rho02/(gamma02-1); S(3,2) = rho01/(gamma01-1); S(3,3) = rho01/(gamma01-1);
#elif GAMNCON									                               
  S(3,0) = 1.0/(gamma4-1);  S(3,1) = 1.0/(gamma02-1);   S(3,2) = 1.0/(gamma01-1);   S(3,3) = 1.0/(gamma01-1);  
#endif
#elif TWOD
  if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);
  scalar v4 = 0;            scalar v02 = 0;             scalar v01 = 0;             scalar v03 = 0;
  Et4  = 1.0/(gamma4-1.0)*p4 + 0.5*rho4*(u4*u4+v4*v4);
  Et02  = 1.0/(gamma02-1.0)*p02 + 0.5*rho02*(u02*u02+v02*v02);
  Et01  = 1.0/(gamma01-1.0)*p01 + 0.5*rho01*(u01*u01+v01*v01);
  Et03  = 1.0/(gamma03-1.0)*p03 + 0.5*rho03*(u03*u03+v03*v03);
  S(0,0) = rho4;            S(0,1) = rho02;             S(0,2) = rho01;             S(0,3) = rho03;		    
  S(1,0) = rho4*u4;         S(1,1) = rho02*u02;         S(1,2) = rho01*u01;	    S(1,3) = rho03*u03;	    
  S(2,0) = rho4*v4;         S(2,1) = rho02*v02;         S(2,2) = rho01*v01;	    S(2,3) = rho03*v03;	    
  S(3,0) = Et4;             S(3,1) = Et02;              S(3,2) = Et01; 		    S(3,3) = Et03; 		    
#ifdef GAMCONS									                                
  S(4,0) = rho4/(gamma4-1); S(4,1) = rho02/(gamma02-1); S(4,2) = rho01/(gamma01-1); S(4,3) = rho03/(gamma03-1); 
#elif GAMNCON									                                
  S(4,0) = 1.0/(gamma4-1);  S(4,1) = 1.0/(gamma02-1);   S(4,2) = 1.0/(gamma01-1);   S(4,3) = 1.0/(gamma03-1);   
#endif
#endif

  int ind = 0; // index for the state we use
  for(int e = 0; e < N_E; e++){
    scalar x = XYZNodes(0,e*D+0);
    for(int i = 0; i < N_s; i++){
      for(int b = 0; b < N_S; b++) if ((xS(b,0) <= x) && (x < xS(b+1,0)))  ind = b;
      for(int k = 0; k < N_F; k++) U(i,e*N_F+k) = S(k,ind);
    }
  }
}

void init_dg_sinephi_passive(const int N_s, const int N_E, const int N_F, const int D, scalar &gamma, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  scalar rho     = 1.0;
  scalar u       = 1.0;
  gamma          = 1.4;
  scalar phi     = 0.0;
  scalar sinephi = 0.0; // sine perturbation on phi
  scalar sinerho = 0.0; // sine perturbation on rho
  scalar Aphi    = 0.2; // amplitude of the perturbation
  scalar Arho    = 0.2; // amplitude of the perturbation
  scalar p       = 1.0;
  scalar Et      = 1.0/(gamma-1.0)*p + 0.5*rho*u*u;

  int N_S = 1; // number of states
  fullMatrix<scalar> xS(N_S+1,1); // state boundaries
  xS(0,0) = -1000; xS(1,0) = 1000;
  fullMatrix<scalar> S(N_F,N_S); // holds the states
  
  int ind = 0; // index for the state we use
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      for(int b = 0; b < N_S; b++) if ((xS(b,0) <= x) && (x < xS(b+1,0)))  ind = b;
      if ((x<-1)||(x>1)) {sinephi = 0; sinerho=0;}
      else               {sinephi = Aphi*sin(M_PI*x); sinerho = Arho*sin(4*M_PI*x);}
#ifdef ONED
      if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);
      S(0,0) = (rho+sinerho);
      S(1,0) = (rho+sinerho)*u;
      S(2,0) = 1.0/(gamma-1.0)*p + 0.5*(rho+sinerho)*u*u;
      S(3,0) = (rho+sinerho)*(phi+sinephi);
      S(4,0) = phi+sinephi;
#elif TWOD
      if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 6.\n",N_F);
      scalar v = 0;            
      S(0,0) = (rho+sinerho);
      S(1,0) = (rho+sinerho)*u;
      S(2,0) = (rho+sinerho)*v;
      S(3,0) = 1.0/(gamma-1.0)*p + 0.5*(rho+sinerho)*(u*u+v*v);
      S(4,0) = (rho+sinerho)*(phi+sinephi);
      S(5,0) = phi+sinephi;
#endif
      for(int k = 0; k < N_F; k++) U(i,e*N_F+k) = S(k,ind);
    }
  }
}


void init_dg_sodmono_passive(const int N_s, const int N_E, const int N_F, const int D, scalar &gamma, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  gamma = 1.4;
  // Left state
  scalar rhoL = 1;
  scalar uL   = 0;
  scalar pL   = 1.0;
  scalar EtL  = 1.0/(gamma-1.0)*pL + 0.5*rhoL*uL*uL;
    
  // Right state
  scalar rhoR = 0.125;
  scalar uR   = 0;
  scalar pR   = 0.1;
  scalar EtR  = 1.0/(gamma-1.0)*pR + 0.5*rhoR*uR*uR;
  
  int N_S = 2; // number of states
  fullMatrix<scalar> xS(N_S+1,1); // state boundaries
  xS(0,0) = -1000; xS(1,0) = 1E-8; xS(2,0) = 1000;
  fullMatrix<scalar> S(N_F,N_S); // holds the states

#ifdef ONED
  if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);
  S(0,0) = rhoL;    S(0,0) = rhoR;   
  S(1,0) = rhoL*uL; S(1,0) = rhoR*uR;
  S(2,0) = EtL;	    S(2,0) = EtR;	   
  S(3,0) = 0.0;	    S(3,0) = 0.0;	   
  S(4,0) = 0.0;     S(4,0) = 0.0;    
#elif TWOD
  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 6.\n",N_F);
  scalar vL = 0;    scalar vR = 0;       
  EtL  = 1.0/(gamma-1.0)*pL + 0.5*rhoL*(uL*uL+vL*vL);
  EtR  = 1.0/(gamma-1.0)*pR + 0.5*rhoR*(uR*uR+vR*vR);
  S(0,0) = rhoL;    S(0,0) = rhoR;   
  S(1,0) = rhoL*uL; S(0,0) = rhoR*uR;   
  S(2,0) = rhoL*vL; S(1,0) = rhoR*vR;
  S(3,0) = EtL;	    S(2,0) = EtR;	   
  S(4,0) = 0.0;	    S(3,0) = 0.0;	   
  S(5,0) = 0.0;	    S(4,0) = 0.0;    
#endif

  int ind = 0; // index for the state we use
  for(int e = 0; e < N_E; e++){
    scalar x = XYZNodes(0,e*D+0);
    for(int i = 0; i < N_s; i++){
      for(int b = 0; b < N_S; b++) if ((xS(b,0) <= x) && (x < xS(b+1,0)))  ind = b;
      for(int k = 0; k < N_F; k++) U(i,e*N_F+k) = S(k,ind);
    }
  }
}

void init_dg_euler1D_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions (see CFD2 project 1)

  scalar a = 0;
  scalar u = 0;
  scalar rho = 0;

  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if (x<=-1.5){
	u = -2.0/gamma;
	a = (1-gamma)/gamma + 1;
	rho = gamma*pow(a,2.0/(gamma-1));
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)0;
	U(i,e*N_F+3) = (scalar)0;
	U(i,e*N_F+4) = (scalar)0;
	U(i,e*N_F+5) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
      }
      else if ((x>-1.5)&&(x<-0.5)){
	u = -1.0/gamma*(1-tanh((x+1)/(0.25-(x+1)*(x+1))));
	a = u*(gamma-1)/2.0 + 1.0;
	rho = gamma*pow(a,2.0/(gamma-1));;
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)0;
	U(i,e*N_F+3) = (scalar)0;
	U(i,e*N_F+4) = (scalar)0;
	U(i,e*N_F+5) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
      }
      else if ((x>=-0.5)&&(x<=0.5)){
	u = 0;
	a = 0;
	rho = gamma;
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)0;
	U(i,e*N_F+3) = (scalar)0;
	U(i,e*N_F+4) = (scalar)0;
	U(i,e*N_F+5) = (scalar)1.0/(gamma-1);
      }
      else if ((x>0.5)&&(x<1.5)){
	u = 1.0/gamma*(1+tanh((x-1)/(0.25-(x-1)*(x-1))));;
	a = 1 - (gamma-1)/2.0*u;
	rho = gamma*pow(a,2.0/(gamma-1));
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)0;
	U(i,e*N_F+3) = (scalar)0;
	U(i,e*N_F+4) = (scalar)0;
	U(i,e*N_F+5) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
      }
      else if (x>=1.5){
	u = 2.0/gamma;
	a = 1 - (gamma-1)/gamma;
	rho = gamma*pow(a,2.0/(gamma-1));
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)0;
	U(i,e*N_F+3) = (scalar)0;
	U(i,e*N_F+4) = (scalar)0;
	U(i,e*N_F+5) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
      }
    }
  }
}

void init_dg_euler2D_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions (see CFD2 project 2)

  scalar q = 0;
  scalar a = 0; // speed of sound
  scalar rho = 0.0;
  scalar u = 0.0;
  scalar v = 0.0;


  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      scalar y = XYZNodes(i,e*D+1);
      scalar rad = sqrt(x*x+y*y);

      if      ((rad>=0.0)&&(rad<0.5)) q = 0.0;
      else if ((rad>=0.5)&&(rad<1.5)) q = 1.0/gamma * (1.0 + tanh((rad-1.0)/(0.25 - (rad-1.0)*(rad-1.0)))) ;
      else                            q = 2.0/gamma;

      a = 1.0 - q*(gamma-1)/2.0;
      rho = gamma*pow(a,2.0/(gamma-1));
      u = q*x/rad;
      v = q*y/rad;
      U(i,e*N_F+0) = (scalar)rho;
      U(i,e*N_F+1) = (scalar)rho*u;
      U(i,e*N_F+2) = (scalar)rho*v;	
      U(i,e*N_F+3) = (scalar)0;
      U(i,e*N_F+4) = (scalar)0;
      U(i,e*N_F+5) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*(rho*u*rho*u + rho*v*rho*v)/rho;
    }
  }
}

void init_dg_sodtube_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions
  // U = (  rho, rho ux, rho uy, rho uz,   Bx, By, Bz,    E,   ee)
  //   = (    1,      0,      0,      0,   0,  0,  0, 1.78,  0.5)  for (x<0)
  //   = (0.125,      0,      0,      0,   0,  0,  0, 0.88, 0.05)  for (x>=0)
  // gamma = 1.4

  // Left state
  scalar rhoL = 1;
  scalar uL   = 0;
  scalar vL   = 0;
  scalar BxL  = 0;
  scalar ByL  = 0;
  scalar pL   = 1.0;
  scalar EtL  = 1.0/(gamma-1.0)*pL + 0.5*(rhoL*(uL*uL + vL*vL) + (BxL*BxL + ByL*ByL));

  // Right state
  scalar rhoR = 0.125;
  scalar uR   = 0;
  scalar vR   = 0;
  scalar BxR  = 0;
  scalar ByR  = 0;
  scalar pR   = 0.1;
  scalar EtR  = 1.0/(gamma-1.0)*pR + 0.5*(rhoR*(uR*uR + vR*vR) + (BxR*BxR + ByR*ByR));
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if (x<0){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = rhoL*vL;
	U(i,e*N_F+3) = BxL ;
	U(i,e*N_F+4) = ByL ;
	U(i,e*N_F+5) = EtL ;
      }
      else if (x>=0){
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = rhoR*vR;
	U(i,e*N_F+3) = BxR ;
	U(i,e*N_F+4) = ByR ;
	U(i,e*N_F+5) = EtR ; 
      }
    }
  }
}


void init_dg_explode_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions (see Toro p. 587)

  // Left state
  scalar rhoL = 1;
  scalar uL   = 0;
  scalar vL   = 0;
  scalar BxL  = 0;
  scalar ByL  = 0;
  scalar pL   = 1.0;
  scalar EtL  = 1.0/(gamma-1.0)*pL + 0.5*(rhoL*(uL*uL + vL*vL) + (BxL*BxL + ByL*ByL));
 
  // Right state
  scalar rhoR = 0.125;
  scalar uR   = 0;
  scalar vR   = 0;
  scalar BxR  = 0;
  scalar ByR  = 0;
  scalar pR   = 0.1;
  scalar EtR  = 1.0/(gamma-1.0)*pR + 0.5*(rhoR*(uR*uR + vR*vR) + (BxR*BxR + ByR*ByR));

  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      scalar y = XYZNodes(i,e*D+1);
      scalar rad = sqrt(x*x+y*y);

      if (rad<0.4){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = rhoL*vL;
	U(i,e*N_F+3) = BxL ;
	U(i,e*N_F+4) = ByL ;
	U(i,e*N_F+5) = EtL ;
      }
      else{
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = rhoR*vR;
	U(i,e*N_F+3) = BxR ;
	U(i,e*N_F+4) = ByR ;
	U(i,e*N_F+5) = EtR ; 
      }
    }
  }
}




void init_dg_brio_wu_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions
  // U = (  rho, rho ux, rho uy, rho uz,   Bx, By, Bz,    E,   ee)
  //   = (    1,      0,      0,      0, 0.75,  1,  0, 1.78,  0.5)  for (x<0)
  //   = (0.125,      0,      0,      0, 0.75, -1,  0, 0.88, 0.05)  for (x>=0)
  // have gamma = 2

  gamma = 2;

  // Left state
  scalar rhoL = 1;
  scalar uL   = 0;
  scalar vL   = 0;
  scalar BxL  = 0.75;
  scalar ByL  = 1;
  scalar pL   = 1.0;
  scalar EtL  = 1.0/(gamma-1.0)*pL + 0.5*(rhoL*(uL*uL + vL*vL) + (BxL*BxL + ByL*ByL));

  // Right state
  scalar rhoR = 0.125;
  scalar uR   = 0;
  scalar vR   = 0;
  scalar BxR  = 0.75;
  scalar ByR  = -1;
  scalar pR   = 0.1;
  scalar EtR  = 1.0/(gamma-1.0)*pR + 0.5*(rhoR*(uR*uR + vR*vR) + (BxR*BxR + ByR*ByR));

  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if (x<0){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = rhoL*vL;
	U(i,e*N_F+3) = BxL ;
	U(i,e*N_F+4) = ByL ;
	U(i,e*N_F+5) = EtL ;
      }
      else if (x>=0){
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = rhoR*vR;
	U(i,e*N_F+3) = BxR ;
	U(i,e*N_F+4) = ByR ;
	U(i,e*N_F+5) = EtR ; 
      }
    }
  }
}

void init_dg_alfvenw_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions paper by Falle
  // http://www.astro.uni-bonn.de/~jmackey/jmac/node7.html
  // have gamma = 5/3

  gamma = 5.0/3.0;

  // Left state
  scalar rhoL = 1;
  scalar uL   = 0;
  scalar vL   = 1;
  scalar BxL  = 1;
  scalar ByL  = 1;
  scalar pL   = 1.0;
  scalar EtL  = 1.0/(gamma-1.0)*pL + 0.5*(rhoL*(uL*uL + vL*vL) + (BxL*BxL + ByL*ByL));

  // Right state
  scalar rhoR = 1;
  scalar uR   = 0;
  scalar vR   = 1;
  scalar BxR  = 1;
  scalar ByR  = 1;
  scalar pR   = 1.0;
  scalar EtR  = 1.0/(gamma-1.0)*pR + 0.5*(rhoR*(uR*uR + vR*vR) + (BxR*BxR + ByR*ByR));


  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if (x<0){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = rhoL*vL;
	U(i,e*N_F+3) = BxL ;
	U(i,e*N_F+4) = ByL ;
	U(i,e*N_F+5) = EtL ;
      }
      else if (x>=0){
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = rhoR*vR;
	U(i,e*N_F+3) = BxR ;
	U(i,e*N_F+4) = ByR ;
	U(i,e*N_F+5) = EtR ; 
      }
    }
  }
}

void init_dg_fastshk_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions paper by Falle
  // http://www.astro.uni-bonn.de/~jmackey/jmac/node7.html
  // have gamma = 5/3

  gamma = 5.0/3.0;

  // Left state
  scalar rhoL = 3;
  scalar uL   = -0.732;
  scalar vL   = -1.333;
  scalar BxL  = 3;
  scalar ByL  = 2.309;
  scalar pL   = 16.33;
  scalar EtL  = 1.0/(gamma-1.0)*pL + 0.5*(rhoL*(uL*uL + vL*vL) + (BxL*BxL + ByL*ByL));

  // Right state
  scalar rhoR = 1;
  scalar uR   = -4.196;
  scalar vR   = 0;
  scalar BxR  = 3;
  scalar ByR  = 0;
  scalar pR   = 1.0;
  scalar EtR  = 1.0/(gamma-1.0)*pR + 0.5*(rhoR*(uR*uR + vR*vR) + (BxR*BxR + ByR*ByR));

  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if (x<0){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = rhoL*vL;
	U(i,e*N_F+3) = BxL ;
	U(i,e*N_F+4) = ByL ;
	U(i,e*N_F+5) = EtL ;
      }
      else if (x>=0){
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = rhoR*vR;
	U(i,e*N_F+3) = BxR ;
	U(i,e*N_F+4) = ByR ;
	U(i,e*N_F+5) = EtR ; 
      }
    }
  }
}

void init_dg_slowshk_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions paper by Falle
  // http://www.astro.uni-bonn.de/~jmackey/jmac/node7.html
  // have gamma = 5/3

  gamma = 5.0/3.0;

  // Left state
  scalar rhoL = 1.368;
  scalar uL   = 0.269;
  scalar vL   = 1;
  scalar BxL  = 1;
  scalar ByL  = 0;
  scalar pL   = 1.769;
  scalar EtL  = 1.0/(gamma-1.0)*pL + 0.5*(rhoL*(uL*uL + vL*vL) + (BxL*BxL + ByL*ByL));

  // Right state
  scalar rhoR = 1;
  scalar uR   = 0;
  scalar vR   = 0;
  scalar BxR  = 1;
  scalar ByR  = 1;
  scalar pR   = 1.0;
  scalar EtR  = 1.0/(gamma-1.0)*pR + 0.5*(rhoR*(uR*uR + vR*vR) + (BxR*BxR + ByR*ByR));

  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if (x<0){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = rhoL*vL;
	U(i,e*N_F+3) = BxL ;
	U(i,e*N_F+4) = ByL ;
	U(i,e*N_F+5) = EtL ;
      }
      else if (x>=0){
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = rhoR*vR;
	U(i,e*N_F+3) = BxR ;
	U(i,e*N_F+4) = ByR ;
	U(i,e*N_F+5) = EtR ; 
      }
    }
  }
}

void init_dg_fastrar_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions paper by Falle
  // http://www.astro.uni-bonn.de/~jmackey/jmac/node7.html
  // have gamma = 5/3

  gamma = 5.0/3.0;

  // Left state
  scalar rhoL = 1;
  scalar uL   = 0;
  scalar vL   = 0;
  scalar BxL  = 1;
  scalar ByL  = 3;
  scalar pL   = 2;
  scalar EtL  = 1.0/(gamma-1.0)*pL + 0.5*(rhoL*(uL*uL + vL*vL) + (BxL*BxL + ByL*ByL));

  // Right state
  scalar rhoR = 0.2641;
  scalar uR   = 3.6;
  scalar vR   = -2.551;
  scalar BxR  = 1;
  scalar ByR  = 0;
  scalar pR   = 0.2175;
  scalar EtR  = 1.0/(gamma-1.0)*pR + 0.5*(rhoR*(uR*uR + vR*vR) + (BxR*BxR + ByR*ByR));
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if (x<0){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = rhoL*vL;
	U(i,e*N_F+3) = BxL ;
	U(i,e*N_F+4) = ByL ;
	U(i,e*N_F+5) = EtL ;
      }
      else if (x>=0){
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = rhoR*vR;
	U(i,e*N_F+3) = BxR ;
	U(i,e*N_F+4) = ByR ;
	U(i,e*N_F+5) = EtR ; 
      }
    }
  }
}

void init_dg_slowrar_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

  if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions paper by Falle
  // http://www.astro.uni-bonn.de/~jmackey/jmac/node7.html
  // have gamma = 5/3

  gamma = 5.0/3.0;

  // Left state
  scalar rhoL = 1;
  scalar uL   = 0;
  scalar vL   = 0;
  scalar BxL  = 1;
  scalar ByL  = 0;
  scalar pL   = 2;
  scalar EtL  = 1.0/(gamma-1.0)*pL + 0.5*(rhoL*(uL*uL + vL*vL) + (BxL*BxL + ByL*ByL));

  // Right state
  scalar rhoR = 0.2;
  scalar uR   = 1.186;
  scalar vR   = 2.967;
  scalar BxR  = 1;
  scalar ByR  = 1.6405;
  scalar pR   = 0.1368;
  scalar EtR  = 1.0/(gamma-1.0)*pR + 0.5*(rhoR*(uR*uR + vR*vR) + (BxR*BxR + ByR*ByR));
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if (x<0){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = rhoL*vL;
	U(i,e*N_F+3) = BxL ;
	U(i,e*N_F+4) = ByL ;
	U(i,e*N_F+5) = EtL ;
      }
      else if (x>=0){
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = rhoR*vR;
	U(i,e*N_F+3) = BxR ;
	U(i,e*N_F+4) = ByR ;
	U(i,e*N_F+5) = EtR ; 
      }
    }
  }
}


void init_dg_ovortex_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

  // http://www.astro.princeton.edu/~jstone/Athena/tests/orszag-tang/pagesource.html
  // http://www.astro.virginia.edu/VITA/ATHENA/ot.html
  // http://flash.uchicago.edu/site/flashcode/user_support/flash4_ug/node33.html#SECTION08122000000000000000

  
  scalar rho = 25.0/(36.0*M_PI);
  scalar p = 5.0/(12*M_PI);
  scalar B0= 1/sqrt(4*M_PI);
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      scalar y = XYZNodes(i,e*D+1);
      scalar u = -rho*sin(2*M_PI*y);
      scalar v = -rho*sin(2*M_PI*x);
      scalar Bx= -B0*sin(2*M_PI*y);
      scalar By=  B0*sin(4*M_PI*x);
      scalar Et= 1.0/(gamma-1.0)*p + 0.5*(rho*(u*u + v*v) + (Bx*Bx + By*By));
      U(i,e*N_F+0) = rho;
      U(i,e*N_F+1) = u;
      U(i,e*N_F+2) = v;
      U(i,e*N_F+3) = Bx;
      U(i,e*N_F+4) = By;
      U(i,e*N_F+5) = Et;
    }
  }
}

void init_dg_mhdroto_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

  // http://flash.uchicago.edu/site/flashcode/user_support/flash4_ug/node33.html#SECTION08122000000000000000
  
  scalar r0 = 0.1;
  scalar r1 = 0.115;
  scalar r  = 0.0;
  scalar fr = 0.0;
  scalar p = 1;
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      scalar y = XYZNodes(i,e*D+1);
      r = sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
      fr = (r1-r)/(r1-r0);
      if (r<=r0){      
	scalar rho = 10;
	scalar u = -fr*(y-0.5)/r0;	 
	scalar v =  fr*(x-0.5)/r0;	 
	scalar Bx= 5.0/sqrt(4*M_PI); 
	scalar By= 0;		 
	scalar Et= 1.0/(gamma-1.0)*p + 0.5*(rho*(u*u + v*v) + (Bx*Bx + By*By));
	U(i,e*N_F+0) = rho;
	U(i,e*N_F+1) = u;
	U(i,e*N_F+2) = v;
	U(i,e*N_F+3) = Bx;
	U(i,e*N_F+4) = By;
	U(i,e*N_F+5) = Et;
      }
      else if ((r>r0)&&(r<r1)){
	scalar rho= 1.0+9.0*fr;	  
	scalar u =  -fr*(y-0.5)/r;	  
	scalar v =   fr*(x-0.5)/r;	  
	scalar Bx=  5.0/sqrt(4*M_PI);  
	scalar By=  0;		  
	scalar Et= 1.0/(gamma-1.0)*p + 0.5*(rho*(u*u + v*v) + (Bx*Bx + By*By));
	U(i,e*N_F+0) = rho;
	U(i,e*N_F+1) = u;
	U(i,e*N_F+2) = v;
	U(i,e*N_F+3) = Bx;
	U(i,e*N_F+4) = By;
	U(i,e*N_F+5) = Et;
      }
      else if (r>=r1){
	scalar rho=1;		  
	scalar u = 0;		  
	scalar v = 0;		  
	scalar Bx= 5.0/sqrt(4*M_PI);  
	scalar By= 0;		  
	scalar Et= 1.0/(gamma-1.0)*p + 0.5*(rho*(u*u + v*v) + (Bx*Bx + By*By));
	U(i,e*N_F+0) = rho;
	U(i,e*N_F+1) = u;
	U(i,e*N_F+2) = v;
	U(i,e*N_F+3) = Bx;
	U(i,e*N_F+4) = By;
	U(i,e*N_F+5) = Et;
      }
    }
  }
}

