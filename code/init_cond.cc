#include <init_cond.h>
#include <stdlib.h>

void init_dg_shallow(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

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


void init_dg_simplew_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

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

void buildLRstates_multifluid(scalar rhoL, scalar uL, scalar EtL, scalar gammaL, scalar rhoR, scalar uR, scalar EtR, scalar gammaR, const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

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

void init_dg_sodtube_multifluid(const int N_s, const int N_E,const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

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
  
  buildLRstates_multifluid(rhoL, uL, EtL, gammaL, rhoR, uR, EtR, gammaR, N_s, N_E, XYZNodes, U);

}

void init_dg_sodmono_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

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
  
  buildLRstates_multifluid(rhoL, uL, EtL, gammaL, rhoR, uR, EtR, gammaR, N_s, N_E, XYZNodes, U);
}


void init_dg_contact_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

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

  buildLRstates_multifluid(rhoL, uL, EtL, gammaL, rhoR, uR, EtR, gammaR, N_s, N_E, XYZNodes, U);

}

void init_dg_rhotact_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){
  
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

  buildLRstates_multifluid(rhoL, uL, EtL, gammaL, rhoR, uR, EtR, gammaR, N_s, N_E, XYZNodes, U);
}


void init_dg_matfrnt_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

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

  scalar GL, GR;
#ifdef GAMCONS
  GL = rhoL/(gammaL-1.0);   GR = rhoR/(gammaR-1.0);
#elif GAMNCON
  GL = 1.0/(gammaL-1.0);    GR = 1.0/(gammaR-1.0);
#endif

  scalar xc=0;
  for(int e = 0; e < N_E; e++){
    xc = XYZCen(e,0);
    for(int i = 0; i < N_s; i++){
      if ((-0.5<xc)&&(xc<0.5)){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = EtL;
	U(i,e*N_F+3) = GL;
      }
      else {
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = EtR;
	U(i,e*N_F+3) = GR;
      }
    }
  }
}

void init_dg_sinegam_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

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
      sinerho = Arho*sin(4*M_PI*x);
      sinegam = Agam*sin(2*M_PI*x);
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
      scalar y = XYZNodes(i,e*D+1);
      sinerho = Arho*sin(2.0*M_PI*x)*sin(2.0*M_PI*y);
      sinegam = Agam*sin(2.0*M_PI*x)*sin(2.0*M_PI*y);
      scalar v = 1;
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
             
void init_dg_expogam_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

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

void init_dg_shckint_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

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

void init_dg_multint_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

#ifdef TWOD
  printf("multint problem can only be run in 1D. Exiting");
  exit(1);
#endif

  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);

  // Initialize
  scalar xshck = 0.05; // initial shock location
  scalar Lx = 0.089*2.0/3.0; // wavelength (or width of tube)
  scalar K = 10; 
  scalar h = K*Lx;
  scalar xinterface1 = 0; // first interface location
  scalar xinterface2 =-h; // second interface location
  scalar delta=0.005;               // The diffusion layer thickness

  // Velocities/pressures in all materials
  scalar ucoord = 0; // coordinate shift to the right
  scalar u = 0.0+ucoord;
  scalar p = 1e5;

  // Convention for material order:
  // mat1 (with shock) | mat 2 | mat 3
  // pre-shock density (material 1)
  // The shock is initialized in here
  scalar rho01   = 1.351;
  scalar gamma01 = 1.276;
  scalar alpha01 = 1/(gamma01-1);
  scalar c01     = sqrt(gamma01*p/rho01); // sound speed
  scalar M1      = 34.76; // molecular weight

  // pre-shock density (material 2)
  scalar rho02   = 5.494;
  scalar gamma02 = 1.093;
  scalar alpha02 = 1/(gamma02-1);
  scalar M2      = 146.05;

  // pre-shock density (material 3)
  scalar rho03   = 0.1785; // 10
  scalar gamma03 = 5.0/3.0;
  scalar alpha03 = 1/(gamma03-1);
  scalar M3      = 300; //4
  
  // Post-shock state (material 1) (see p 101 Toro)
  scalar Ms     = 1.21;   // Shock Mach number
  scalar u4     =-Ms*c01*(2*(Ms*Ms-1))/(gamma01+1)/(Ms*Ms)+ucoord; // shock is moving downwards
  scalar p4     = p*(1+2*gamma01/(gamma01+1)*(Ms*Ms-1));
  scalar rho4   = rho01*(gamma01+1)*Ms*Ms/(2+(gamma01-1)*Ms*Ms);
  scalar gamma4 = gamma01;
  scalar Et4    = p4/(gamma4-1.0) + 0.5*rho4*u4*u4;

  for(int e = 0; e < N_E; e++){
    scalar xc = XYZCen(e,0);
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);

      if(xc >= (xshck-1e-6)){ // post-shock region
	U(i,e*N_F+0) = rho4;
	U(i,e*N_F+1) = rho4*u4;
	U(i,e*N_F+2) = Et4 ;
#ifdef GAMCONS
	U(i,e*N_F+3) = rho4/(gamma4-1);
#elif GAMNCON
	U(i,e*N_F+3) = 1.0/(gamma4-1);
#endif
      }
      else{
  	// horizontal distance from interface
	scalar d1 = (delta-x+xinterface1)/(2*delta);
	scalar d2 = (delta-x+xinterface2)/(2*delta);
	
	// Calculate volume fractions
	scalar vol1=0;
	scalar vol2=0;
	if      (d1<=0)            { vol1 = 1; vol2 = 0;}
	else if ((0<d1)&&(d1<1))   { vol1 = exp(log(1e-16)*pow(fabs(d1),8)); vol2 = 1-vol1;}
	else{
	  if      (d2<=0)          { vol1 = 0; vol2 = 1;}
	  else if ((0<d2)&&(d2<1)) { vol1 = 0; vol2 = exp(log(1e-16)*pow(fabs(d2),8));}
	  else                     { vol1 = 0; vol2 = 0;}  
	}
	
	scalar j1  = vol1;
	scalar j2  = vol2;
	scalar j3  = 1-vol1-vol2;
	scalar rho = j1*rho01+j2*rho02+j3*rho03;
	scalar Y1  = j1*rho01/rho;       // mass fraction
	scalar Y2  = j2*rho02/rho;
	scalar Y3  = j3*rho03/rho;
	scalar M   = 1/(Y1/M1+Y2/M2+Y3/M3);                    // total molecular weight

	scalar alpha = Y1*alpha01*M/M1+Y2*alpha02*M/M2+Y3*alpha03*M/M3;
	scalar gamma = 1+1.0/alpha;
	
	U(i,e*N_F+0) = rho;
	U(i,e*N_F+1) = rho*u;
	U(i,e*N_F+2) = p/(gamma-1)+ 0.5*rho*u*u;
#ifdef GAMCONS
	U(i,e*N_F+3) = rho/(gamma-1);
#elif GAMNCON
	U(i,e*N_F+3) = 1.0/(gamma-1);
#endif
      }
    }
  }
}

void init_dg_blast1d_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

#ifdef TWOD
  printf("blast1d problem can only be run in 1D. Exiting");
  exit(1);
#endif

  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);

  scalar gamma = 5./3.;
  scalar alpha = 2./3.; // = 2/3 for planar, 1/2 for cyl, 2/5 for sph.
  scalar Q = 0.66927; // for alpha = 2/3 and gamma = 5/3

  // Initialize by setting the explosion energy
  scalar patm = 1e5;
  scalar u0  = 0;    // velocity of unshocked material

  //scalar rho0 = 100; // density of unshocked material (100 kg/m^3 = 0.1 g/cc)
  //scalar p0 = 1.98e11; //1.6e11  // pressure at shock in Pa
  //scalar t0 = 35*1e-9; // time = 25ns
  //scalar R0 = sqrt(0.5*(gamma+1)*p0/rho0)/(alpha*pow(t0,alpha-1));
  //scalar Ex = rho0*pow(Q,3)*pow(R0,3); // explosion energy
  scalar rho0 = 50;
  scalar Ex = 1.67e8;
  printf("Explosion energy=%e\n",Ex);
  scalar blastpos = 0.0;
  scalar Dxx = 0.00005;
  
  for(int e = 0; e < N_E; e++){
    scalar xc = XYZCen(e,0);
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);

      if(xc<blastpos){
	U(i,e*N_F+0) = rho0;
	U(i,e*N_F+1) = rho0*u0;
	U(i,e*N_F+2) = Ex/Dxx;
#ifdef GAMCONS
	U(i,e*N_F+3) = rho0/(gamma-1);
#elif GAMNCON
	U(i,e*N_F+3) = 1.0/(gamma-1);
#endif
      }
      else{
	U(i,e*N_F+0) = rho0;
	U(i,e*N_F+1) = rho0*u0;
	U(i,e*N_F+2) = patm/(gamma-1.0) + 0.5*rho0*u0*u0;
#ifdef GAMCONS
	U(i,e*N_F+3) = rho0/(gamma-1);
#elif GAMNCON
	U(i,e*N_F+3) = 1.0/(gamma-1);
#endif
      }
    }
  }

//   // Initialize using interpolation method
  // // Load/open the table file
  // std::string tablefile = "xwup.txt";

  // // Read the data from file and fill a matrix xwup
  // fullMatrix<scalar> XWUP;
  // scalar gamma;
  // scalar alpha;
  // scalar Q;
  // readTable(tablefile.c_str(),XWUP,gamma,alpha,Q);

//   scalar patm = 1e5;
//   scalar p0 = 1.6e11;  // pressure at shock in Pa
//   scalar t0 = 25*1e-9; // time = 25ns
//   scalar rho0 = 100; // density of unshocked material (100 kg/m^3 = 0.1 g/cc)
//   scalar u0  = 0;    // velocity of unshocked material
//   scalar R0 = sqrt(0.5*(gamma+1)*p0/rho0)/(alpha*pow(t0,alpha-1));
  
//   std::vector<std::pair<scalar, scalar> > rho_r;
//   std::vector<std::pair<scalar, scalar> > u_r;
//   std::vector<std::pair<scalar, scalar> > p_r;
//   scalar t = 1e-9;
//   scalar R    = R0*pow(t,alpha);
//   scalar Rdot = alpha*R0*pow(t,alpha-1);

//   for(int k=0; k<XWUP.size1(); k++){

//     scalar r_t = R*XWUP(k,0);
//     scalar rho_t = rho0*XWUP(k,1);
//     scalar u_t = Rdot*XWUP(k,2);
//     scalar p_t = rho0*Rdot*Rdot*XWUP(k,3);
    
//     rho_r.push_back(std::make_pair(r_t,rho_t));
//     u_r.push_back(std::make_pair(r_t,u_t));
//     p_r.push_back(std::make_pair(r_t,p_t));
//   }

  
//   // Blast wave profile
//   for(int e = 0; e < N_E; e++){
//     for(int i = 0; i < N_s; i++){
//       scalar x = XYZNodes(i,e*D+0);

//       scalar rho = interpolate(x,rho_r,rho_r[0].second,rho0);
//       scalar u   = interpolate(x,u_r,0,0);
//       scalar p   = interpolate(x,p_r,p_r[0].second,patm);
      
//       U(i,e*N_F+0) = rho;
//       U(i,e*N_F+1) = rho*u;
//       U(i,e*N_F+2) = p/(gamma-1.0) + 0.5*rho*u*u;
// #ifdef GAMCONS
//       U(i,e*N_F+3) = rho/(gamma-1);
// #elif GAMNCON
//       U(i,e*N_F+3) = 1.0/(gamma-1);
// #endif
//     }
//   }
}

void init_dg_simblst_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

  // Initialize a rarefaction moving towards the left (or downwards)
  scalar A0 = 0.00183;
  scalar Lx = 0.089*2.0/3.0;
  scalar K = 2;
  scalar yinterface =-(K)*Lx;//-0.07; // first interface location
  scalar delta=0.005;   // The diffusion layer thickness
  scalar u0 = 0;
  scalar v0 = 0;
  scalar p0 = 1e5;

  // Parameters to choose
  scalar L = K*Lx;       // length of rarefaction when it reaches the interface
  scalar Ms = 2.85;      // shock mach number
  scalar T = 0.;         // time difference btw shock and rarefaction arrival at the interface
  
  // Top material (shock and rarefaction initialized here)
  scalar rho01   = 1.351;//5.494;//
  scalar gamma01 = 1.093;//1.276;//;
  scalar alpha01 = 1/(gamma01-1);
  scalar M01     = 34.76;//146.05;//
  scalar c01     = sqrt(gamma01*p0/rho01);
  
  // Bottom material (material 2)
  scalar rho02   = 5.494;//1.351;//
  scalar gamma02 = 1.276;//1.093;//
  scalar alpha02 = 1/(gamma02-1);
  scalar M02     = 146.05;//34.76;//  // molecular weight
  scalar c02     = sqrt(gamma02*p0/rho02);

  // Post-shock material
  scalar us     =-Ms*c01; // shock velocity
  scalar u4     = 0;
  scalar v4     =-Ms*c01*(2*(Ms*Ms-1))/(gamma01+1)/(Ms*Ms); // shock is moving downwards
  scalar p4     = p0*(1+2*gamma01/(gamma01+1)*(Ms*Ms-1));
  scalar rho4   = rho01*(gamma01+1)*Ms*Ms/(2+(gamma01-1)*Ms*Ms);
  scalar gamma4 = gamma01;
  scalar Et4    = p4/(gamma4-1.0) + 0.5*rho4*(u4*u4+v4*v4);
  scalar c4     = sqrt(gamma4*p4/rho4);
  
  // contact velocity (velocity behind the rarefaction, positive)
  scalar vc = 0;//2.0*c4/(gamma01-1) * (1 - pow(strength, (gamma01-1)/(2.0*gamma01))) + v4;

  // time at which the rarefaction is of length L
  scalar tiR = -2.0/(gamma01+1) * 1.0/v4 * L;

  // solve for the origin of rarefaction
  scalar yR = yinterface + (c4-v4)*tiR;

  // Time at which the shock gets to the interface
  scalar tiS = tiR - T;

  // Solve for the shock origin
  scalar yS = yinterface - us*tiS;

  // Initialization time t0 = alpha*tiS
  scalar t0 = 0.95*tiS;

  // Position of the shock, head, and tail at t0
  scalar yS0 = us*t0 + yS;
  scalar yH0 = -(c4-v4)*t0 + yR;    
  scalar yT0 = -(c4 - v4 + (gamma01+1)/2.0 * v4) * t0 + yR;
  printf("c4=%f, us=%f, yi=%f, yS=%f, yS0=%f, yH0=%f, yT0=%f, %f, L=%f\n", c4, us, yinterface, yS, yS0, yH0, yT0,yT0-yH0, L);
  
  // reflection coefficient
  //scalar R = (1- rho01*c01/(rho02*c02))/(1+rho01*c01/(rho02*c02));
  scalar vcoord = 0;//-(1-R)*vc; // coordinate shift upwards
  
  scalar xc=0, yc=0, x=0, y=0, u, v, rho, p;
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){

#ifdef ONED
      A0 = 0;
      yc = XYZCen(e,0);
      y  = XYZNodes(i,e*D+0);
#elif TWOD
      xc = XYZCen(e,0);
      x  = XYZNodes(i,e*D+0);
      yc = XYZCen(e,1);
      y  = XYZNodes(i,e*D+1);
#endif

      // vertical distance from interface
      scalar d = ((delta+A0*sin(2*M_PI*x/Lx-M_PI/2))-y+yinterface)/(2*delta);
      
      // Calculate volume fractions
      scalar vol=0;
      if      ((d<1)&&(d>0)) vol = exp(log(1e-16)*pow(fabs(d),8));
      else if (d<=0)         vol = 1;
      else                   vol = 0;
      
      scalar jx  = 1-vol;
      scalar rho = jx*rho02+(1-jx)*rho01;
      scalar jy  = jx*rho02/(jx*rho02+(1-jx)*rho01);      // mass fraction
      scalar jM  = 1/(jy/M02+(1-jy)/M01);                 // total molecular weight
      
      scalar alpha = jy*alpha02*jM/M02+(1-jy)*alpha01*jM/M01;
      scalar gamma = 1+1.0/alpha;
      u=0;v=0;p=p0;

      // Velocity, density and pressure modifications for shock and rarefaction
      if ((yS0 <= yc) && (yc <= yH0)){ // region behind the shock,  in front of the rarefaction
      	u = u4;
      	v = v4;
      	rho = rho4;
      	p   = p4;
      }
      else if ((yH0 < yc) && (yc < yT0)){ // inside rarefaction
      	u = 0;
      	v = 2.0/(gamma+1)*(c4 - v4 + (y-yR)/t0) + v4;
	scalar yp = y - v4*t0;
	scalar vp = 2.0/(gamma+1)*(c4 + (yp-yR)/t0);
      	rho = rho4 * pow(1 - (gamma-1)/2.0 * fabs(vp/c4), 2.0/(gamma-1));
      	p   = p4   * pow(1 - (gamma-1)/2.0 * fabs(vp/c4), 2.0*gamma/(gamma-1));
	//v = v + v4;
      }
      else if (yT0 <= yc){ // behind rarefaction
      	u = 0;
      	v = vc;
	scalar vp = -v4;	  
      	rho = rho4 * pow(1 - (gamma-1)/2.0 * fabs(vp/c4), 2.0/(gamma-1));
      	p   = p4   * pow(1 - (gamma-1)/2.0 * fabs(vp/c4), 2.0*gamma/(gamma-1));
      }
      v = v + vcoord;
      
#ifdef ONED
      U(i,e*N_F+0) = rho;
      U(i,e*N_F+1) = rho*v;
      U(i,e*N_F+2) = p/(gamma-1)+ 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
      U(i,e*N_F+3) = rho/(gamma-1);
#elif GAMNCON
      U(i,e*N_F+3) = 1.0/(gamma-1);
#endif
      // Mass fractions
      U(i,e*N_F+4) = (1-jx)*rho;
      
#elif TWOD
      U(i,e*N_F+0) = rho;
      U(i,e*N_F+1) = rho*u;
      U(i,e*N_F+2) = rho*v;
      U(i,e*N_F+3) = p/(gamma-1)+ 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
      U(i,e*N_F+4) = rho/(gamma-1);
#elif GAMNCON
      U(i,e*N_F+4) = 1.0/(gamma-1);
#endif
      // Mass fractions
      U(i,e*N_F+5) = (1-jx)*rho;
#endif
    }
  }
}

void init_dg_sodcirc_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

#ifdef ONED
  printf("sodcirc problem can only be run in 2D. Exiting");
  exit(1);
#elif TWOD
  
  if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);
  
  // Initial conditions (see Toro p. 587)

  // Left state
  scalar rhoL   = 1;
  scalar uL     = 0;
  scalar vL     = 0;
  scalar pL     = 1.0;
  scalar gammaL = 1.4;
  scalar EtL    = 1.0/(gammaL-1.0)*pL + 0.5*rhoL*(uL*uL + vL*vL);
 
  // Right state
  scalar rhoR   = 0.125;
  scalar uR     = 0;
  scalar vR     = 0;
  scalar pR     = 0.1;
  scalar gammaR = 1.4;
  scalar EtR    = 1.0/(gammaR-1.0)*pR + 0.5*rhoR*(uR*uR + vR*vR);

  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      scalar y = XYZNodes(i,e*D+1);
      scalar rad = sqrt(x*x+y*y);

      if (rad<0.4){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = rhoL*vL;
	U(i,e*N_F+3) = EtL ;
#ifdef GAMCONS
	U(i,e*N_F+4) = rhoL/(gammaL-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gammaL-1);
#endif
      }
      else{
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = rhoR*vR;
	U(i,e*N_F+3) = EtR ;
#ifdef GAMCONS
	U(i,e*N_F+4) = rhoR/(gammaR-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gammaR-1);
#endif
      }
    }
  }
#endif
}

void init_dg_rminstb_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

#ifdef ONED
  printf("rminstb problem can only be run in 2D. Exiting");
  exit(1);
#elif TWOD
  
  if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);

  // Initialize
  scalar A0 = 0.00183;//0.00183;                 // initial amplitude
  scalar yshck = 0.05; // initial shock location
  scalar Lx = 0.089*2.0/3.0;

  scalar vcoord = 0;//72; // coordinate shift upwards
    
  // Velocities/pressures in both materials
  scalar u = 0.0;
  scalar v = 0.0+vcoord;
  scalar p = 1e5;
  
  // pre-shock density (material 1)
  scalar rho01   = 5.494;
  scalar gamma01 = 1.093;
  scalar alpha01 = 1/(gamma01-1);
  scalar M1      = 146.05;         // Molecular weight
  
  // pre-shock density (material 2)
  // The shock is initialized in here
  scalar rho02   = 1.351;
  scalar gamma02 = 1.276;
  scalar alpha02 = 1/(gamma02-1);
  scalar c02     = sqrt(gamma02*p/rho02); // sound speed
  scalar Ma02    = v/c02; //Mach number ahead of shock
  scalar M2      = 34.76;  // molecular weight

  // The diffusion layer thickness
  scalar delta=0.005;
  
  // Post-shock state (material 2) (see p 101 Toro)
  scalar Ms = 1.21;   // Shock Mach number
  scalar u4   = 0;
  scalar v4   =-Ms*c02*(2*(Ms*Ms-1))/(gamma02+1)/(Ms*Ms)+vcoord; // shock is moving downwards
  scalar p4   = p*(1+2*gamma02/(gamma02+1)*(Ms*Ms-1));
  scalar rho4 = rho02*(gamma02+1)*Ms*Ms/(2+(gamma02-1)*Ms*Ms);
  scalar gamma4 = gamma02;
  scalar Et4    = p4/(gamma4-1.0) + 0.5*rho4*(u4*u4+v4*v4);

  for(int e = 0; e < N_E; e++){
    scalar xc = XYZCen(e,0);
    scalar yc = XYZCen(e,1);
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      scalar y = XYZNodes(i,e*D+1);

      if(yc >= yshck){ // post-shock region
	U(i,e*N_F+0) = rho4;
	U(i,e*N_F+1) = rho4*u4;
	U(i,e*N_F+2) = rho4*v4;
	U(i,e*N_F+3) = Et4 ;
#ifdef GAMCONS
	U(i,e*N_F+4) = rho4/(gamma4-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma4-1);
#endif
      }
      else{
	// vertical distance from interface
	scalar d = ((delta+A0*sin(2*M_PI*x/Lx-M_PI/2))-y)/(2*delta);
	
	// Calculate volume fractions
	scalar vol=0;
	if      ((d<1)&&(d>0)) vol = exp(log(1e-16)*pow(fabs(d),8));
	else if (d<=0)         vol = 1;
	else                   vol = 0;

	scalar jx  = 1-vol;
	scalar rho = jx*rho01+(1-jx)*rho02;
	scalar jy  = jx*rho01/(jx*rho01+(1-jx)*rho02);       // mass fraction
	scalar jM  = 1/(jy/M1+(1-jy)/M2);                    // total molecular weight

	scalar alpha = jy*alpha01*jM/M1+(1-jy)*alpha02*jM/M2;
	scalar gamma = 1+1.0/alpha;
	
	//rho = rho02;
	//gamma = gamma02;
	// u = 0;
	// v = 0;
	U(i,e*N_F+0) = rho;
	U(i,e*N_F+1) = rho*u;
	U(i,e*N_F+2) = rho*v;
	U(i,e*N_F+3) = p/(gamma-1)+ 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
	U(i,e*N_F+4) = rho/(gamma-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma-1);
#endif
      }
    }
  }
#endif
}

void init_dg_rmmulti_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){
  
  // Initialize
  scalar A01 = 0.00183;                 // initial amplitude
  scalar A02 = 0.00183;                 // initial amplitude
  scalar yshck = 0.025; // initial shock location
  scalar Lx = 0.089*2.0/3.0;
  scalar K = 1;
  scalar h = K*Lx;
  scalar yinterface1 = 0; // first interface location
  scalar yinterface2 =-h; // second interface location
  scalar delta=0.005;    // The diffusion layer thickness
    
  // Velocities/pressures in all materials
  scalar vcoord = 51.5;//111;//51.5;//134;//72.9; // coordinate shift upwards
  scalar u = 0.0;
  scalar v = 0.0+vcoord;
  scalar p = 1e5;

  // Convention for material order:
  // mat1 (with shock) | mat 2 | mat 3
  // pre-shock density (material 1)
  // The shock is initialized in here
  scalar rho01   = 1.351;//5.494;//1.351;
  scalar gamma01 = 1.276;//1.093;//1.276;
  scalar alpha01 = 1/(gamma01-1);
  scalar c01     = sqrt(gamma01*p/rho01); // sound speed
  scalar M1      = 34.76; // molecular weight

  // pre-shock density (material 2)
  scalar rho02   = 5.494;//1.351;//
  scalar gamma02 = 1.093;//1.276;//
  scalar alpha02 = 1/(gamma02-1);
  scalar M2      = 146.05;//34.76;//

  // pre-shock density (material 3)
  scalar rho03   = 10;//0.1785;//10;//5.494;//10;//
  scalar gamma03 = 5.0/3.0;//1.093;//
  scalar alpha03 = 1/(gamma03-1);
  scalar M3      = 300;//4;//300;//146.05;//300;//

  // Post-shock state (material 1) (see p 101 Toro)
  scalar Ms = 1.21;   // Shock Mach number
  scalar u4   = 0;
  scalar v4   =-Ms*c01*(2*(Ms*Ms-1))/(gamma01+1)/(Ms*Ms)+vcoord; // shock is moving downwards
  scalar p4   = p*(1+2*gamma01/(gamma01+1)*(Ms*Ms-1));
  scalar rho4 = rho01*(gamma01+1)*Ms*Ms/(2+(gamma01-1)*Ms*Ms);
  scalar gamma4 = gamma01;
  scalar Et4    = p4/(gamma4-1.0) + 0.5*rho4*(u4*u4+v4*v4);

  scalar xc=0, yc=0, x=0, y=0;
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){

#ifdef ONED
      A01 = 0; A02 = 0;
      yc = XYZCen(e,0);
      y  = XYZNodes(i,e*D+0);
#elif TWOD
      xc = XYZCen(e,0);
      x  = XYZNodes(i,e*D+0);
      yc = XYZCen(e,1);
      y  = XYZNodes(i,e*D+1);
#endif

      if(yc >= (yshck-1e-6)){ // post-shock region

#ifdef ONED
	U(i,e*N_F+0) = rho4;
	U(i,e*N_F+1) = rho4*v4;
	U(i,e*N_F+2) = Et4 ;
#ifdef GAMCONS
	U(i,e*N_F+3) = rho4/(gamma4-1);
#elif GAMNCON
	U(i,e*N_F+3) = 1.0/(gamma4-1);
#endif
 	// Mass fractions
	U(i,e*N_F+4) = 1*rho4;
	U(i,e*N_F+5) = 0;
#elif TWOD
	U(i,e*N_F+0) = rho4;
	U(i,e*N_F+1) = rho4*u4;
	U(i,e*N_F+2) = rho4*v4;
	U(i,e*N_F+3) = Et4 ;
#ifdef GAMCONS
	U(i,e*N_F+4) = rho4/(gamma4-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma4-1);
#endif
 	// Mass fractions
	U(i,e*N_F+5) = 1*rho4;
	U(i,e*N_F+6) = 0;
#endif 
      }
      else{
	// vertical distance from interface
	scalar d1 = ((delta+A01*sin(2*M_PI*x/Lx-M_PI/2))-y+yinterface1)/(2*delta);
	scalar d2 = ((delta+A02*sin(2*M_PI*x/Lx-M_PI/2))-y+yinterface2)/(2*delta);
	
	// Calculate volume fractions
	scalar vol1=0;
	scalar vol2=0;
	if      (d1<=0)            { vol1 = 1; vol2 = 0;}
	else if ((0<d1)&&(d1<1))   { vol1 = exp(log(1e-16)*pow(fabs(d1),8)); vol2 = 1-vol1;}
	else{
	  if      (d2<=0)          { vol1 = 0; vol2 = 1;}
	  else if ((0<d2)&&(d2<1)) { vol1 = 0; vol2 = exp(log(1e-16)*pow(fabs(d2),8));}
	  else                     { vol1 = 0; vol2 = 0;}  
	}
	
	scalar j1  = vol1;
	scalar j2  = vol2;
	scalar j3  = 1-vol1-vol2;
	scalar rho = j1*rho01+j2*rho02+j3*rho03;
	scalar Y1  = j1*rho01/rho;       // mass fraction
	scalar Y2  = j2*rho02/rho;
	scalar Y3  = j3*rho03/rho;
	scalar M   = 1/(Y1/M1+Y2/M2+Y3/M3);                    // total molecular weight

	scalar alpha = Y1*alpha01*M/M1+Y2*alpha02*M/M2+Y3*alpha03*M/M3;
	scalar gamma = 1+1.0/alpha;

#ifdef ONED
	U(i,e*N_F+0) = rho;
	U(i,e*N_F+1) = rho*v;
	U(i,e*N_F+2) = p/(gamma-1)+ 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
	U(i,e*N_F+3) = rho/(gamma-1);
#elif GAMNCON
	U(i,e*N_F+3) = 1.0/(gamma-1);
#endif
 	// Mass fractions
	U(i,e*N_F+4) = Y1*rho;
	U(i,e*N_F+5) = Y3*rho;
#elif TWOD
	U(i,e*N_F+0) = rho;
	U(i,e*N_F+1) = rho*u;
	U(i,e*N_F+2) = rho*v;
	U(i,e*N_F+3) = p/(gamma-1)+ 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
	U(i,e*N_F+4) = rho/(gamma-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma-1);
#endif
 	// Mass fractions
	U(i,e*N_F+5) = Y1*rho;
	U(i,e*N_F+6) = Y3*rho;
#endif
      }
    }
  }
}


void init_dg_khinstb_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

#ifdef ONED
  printf("khinstb problem can only be run in 2D. Exiting");
  exit(1);
#elif TWOD
  
  if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);

  // Initialize
  // Velocities/pressures in both materials
  scalar u = 0.0;
  scalar v = 0.0;
  scalar p = 1e5;
  scalar shckpos = 0.0015;
  
  // pre-shock density (material 1)
  scalar rho01   = 1400;//5.494;//2800;
  scalar gamma01 = 5.0/3.0;
  scalar alpha01 = 1/(gamma01-1);
  
  // pre-shock density (material 2)
  // The shock is initialized in here
  scalar rho02   = 100;//1.351;//100;
  scalar gamma02 = 5.0/3.0;
  scalar alpha02 = 1/(gamma02-1);
  scalar c02     = sqrt(gamma02*p/rho02); // sound speed

  // Post-shock state (material 2) (see p 101 Toro)
  scalar Ms     = 100;//1.21;//100//   // Shock Mach number
  scalar u4     = 0;
  scalar v4     = -c02*2*(Ms*Ms-1)/((gamma02+1)*Ms); // shock is moving downwards
  scalar p4     = p*(2.0*gamma02*Ms*Ms - (gamma02-1))/(gamma02+1);
  scalar rho4   = rho02*(gamma02+1)*Ms*Ms/(2+(gamma02-1)*Ms*Ms);
  scalar gamma4 = gamma02;
  scalar Et4    = p4/(gamma4-1.0) + 0.5*rho4*(u4*u4+v4*v4);

  printf("rho4=%g, v4=%g, p4=%g, upa=%f\n",rho4,v4,p4,fabs(v4)+sqrt(gamma02*p4/rho4));
  
  for(int e = 0; e < N_E; e++){
    scalar xc = XYZCen(e,0);
    scalar yc = XYZCen(e,1);
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      scalar y = XYZNodes(i,e*D+1);
      if((yc>(shckpos-1e-6))&&(xc>0)){ // shocked region
	U(i,e*N_F+0) = rho4;
	U(i,e*N_F+1) = rho4*u4;
	U(i,e*N_F+2) = rho4*v4;
	U(i,e*N_F+3) = Et4 ;
#ifdef GAMCONS
	U(i,e*N_F+4) = rho4/(gamma4-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma4-1);
#endif
      }
      else if((yc<=shckpos)&&(xc>0)){ // top material, unshocked
	U(i,e*N_F+0) = rho02;
	U(i,e*N_F+1) = rho02*u;
	U(i,e*N_F+2) = rho02*v;
	U(i,e*N_F+3) = p/(gamma02-1) + 0.5*rho02*(u*u+v*v);
#ifdef GAMCONS
	U(i,e*N_F+4) = rho02/(gamma02-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma02-1);
#endif
      }
      else if(xc<=0){// lower material, unshocked
	U(i,e*N_F+0) = rho01;
	U(i,e*N_F+1) = rho01*u;
	U(i,e*N_F+2) = rho01*v;
	U(i,e*N_F+3) = p/(gamma01-1) + 0.5*rho01*(u*u+v*v);
#ifdef GAMCONS
	U(i,e*N_F+4) = rho01/(gamma01-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma01-1);
#endif
      }
    }
  }
#endif
}

void init_dg_khblast_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

#ifdef ONED
  printf("khblast problem can only be run in 2D. Exiting");
  exit(1);
#elif TWOD
  
  if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);

  // Initialize
  // Velocities/pressures in both materials
  scalar u = 0.0;
  scalar v = 0.0;
  scalar p = 1e5;
  scalar blstpos = 0.0;
  scalar gamma = 5./3.;
  scalar alpha = 2./3.; // = 2/3 for planar, 1/2 for cyl, 2/5 for sph.
  scalar Q = 0.66927; // for alpha = 2/3 and gamma = 5/3

  // pre-shock density (material 1)
  scalar rho01   = 1400;
  scalar gamma01 = gamma;
  scalar alpha01 = 1/(gamma01-1);
  
  // pre-shock density (material 2)
  // The blast is initialized in here
  scalar rho02   = 50;
  scalar gamma02 = gamma;
  scalar alpha02 = 1/(gamma02-1);

  // Initialize by setting the explosion energy
  // scalar ps = 1.98e11;  // pressure at shock in Pa
  // scalar t0 = 25*1e-9; // time = 25ns
  // scalar R0 = sqrt(0.5*(gamma02+1)*ps/rho02)/(alpha*pow(t0,alpha-1));
  // scalar Ex  = rho02*pow(Q,3)*pow(R0,3); // explosion energy
  scalar Ex = 1.87e8;
  printf("Explosion energy=%e\n",Ex);
  scalar Dxx = 0.00005; // energy initially deposited in Dxx
  
  for(int e = 0; e < N_E; e++){
    scalar xc = XYZCen(e,0);
    scalar yc = XYZCen(e,1);
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      scalar y = XYZNodes(i,e*D+1);
      if((yc>blstpos)&&(xc>0)){ // blast region
	U(i,e*N_F+0) = rho02;
	U(i,e*N_F+1) = rho02*u;
	U(i,e*N_F+2) = rho02*v;
	U(i,e*N_F+3) = Ex/Dxx;
#ifdef GAMCONS
	U(i,e*N_F+4) = rho02/(gamma02-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma02-1);
#endif
      }
      else if((yc<=blstpos)&&(xc>0)){ // top material, unshocked
	U(i,e*N_F+0) = rho02;
	U(i,e*N_F+1) = rho02*u;
	U(i,e*N_F+2) = rho02*v;
	U(i,e*N_F+3) = p/(gamma02-1) + 0.5*rho02*(u*u+v*v);
#ifdef GAMCONS
	U(i,e*N_F+4) = rho02/(gamma02-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma02-1);
#endif
      }
      else if(xc<=0){// lower material, unshocked
	U(i,e*N_F+0) = rho01;
	U(i,e*N_F+1) = rho01*u;
	U(i,e*N_F+2) = rho01*v;
	U(i,e*N_F+3) = p/(gamma01-1) + 0.5*rho01*(u*u+v*v);
#ifdef GAMCONS
	U(i,e*N_F+4) = rho01/(gamma01-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma01-1);
#endif
      }
    }
  }

#endif
}

void init_dg_khpertu_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

#ifdef ONED
  printf("khpertu problem can only be run in 2D. Exiting");
  exit(1);
#elif TWOD
  
  if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);

  // Initialize
  // Velocities/pressures in both materials
  scalar A0 = 0.00183;
  scalar Lx = 0.089*2.0/3.0;
  scalar xinterface = 0.0;
  scalar delta = 0.005;
  scalar u = 0.0;
  scalar v = 0.0;
  scalar p = 1e5;

  // pre-shock density (material 1)
  scalar rho01   = 5.494;
  scalar gamma01 = 1.093;
  scalar alpha01 = 1/(gamma01-1);
  scalar M1      = 146.05;
  
  // pre-shock density (material 2)
  // The blast is initialized in here
  scalar rho02   = 1.351;
  scalar gamma02 = 1.276;
  scalar alpha02 = 1/(gamma02-1);
  scalar c02     = sqrt(gamma02*p/rho02); // sound speed
  scalar M2      = 34.76; // molecular weight
  
  // Shock stuff
  scalar shckpos = 0.0;
  scalar Ms = 1.21;   // Shock Mach number
  scalar u4   = 0;
  scalar v4   =-Ms*c02*(2*(Ms*Ms-1))/(gamma02+1)/(Ms*Ms); // shock is moving downwards
  scalar p4   = p*(1+2*gamma02/(gamma02+1)*(Ms*Ms-1));
  scalar rho4 = rho02*(gamma02+1)*Ms*Ms/(2+(gamma02-1)*Ms*Ms);
  scalar gamma4 = gamma02;
  scalar Et4    = p4/(gamma4-1.0) + 0.5*rho4*(u4*u4+v4*v4);

  for(int e = 0; e < N_E; e++){
    scalar xc = XYZCen(e,0);
    scalar yc = XYZCen(e,1);
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      scalar y = XYZNodes(i,e*D+1);
      if((yc>shckpos)&&(xc>0)){ // blast region
	U(i,e*N_F+0) = rho4;
	U(i,e*N_F+1) = rho4*u4;
	U(i,e*N_F+2) = rho4*v4;
	U(i,e*N_F+3) = Et4;
#ifdef GAMCONS
	U(i,e*N_F+4) = rho4/(gamma4-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma4-1);
#endif
      }
      else{
	// horizontal distance from interface
	scalar d = ((delta+A0*sin(2*M_PI*y/Lx-M_PI/2))-x+xinterface)/(2*delta);
		// Calculate volume fractions
	scalar vol=0;
	if      ((d<1)&&(d>0)) vol = exp(log(1e-16)*pow(fabs(d),8));
	else if (d<=0)         vol = 1;
	else                   vol = 0;

	scalar jx  = 1-vol;
	scalar rho = jx*rho01+(1-jx)*rho02;
	scalar jy  = jx*rho01/(jx*rho01+(1-jx)*rho02);       // mass fraction
	scalar jM  = 1/(jy/M1+(1-jy)/M2);                    // total molecular weight

	scalar alpha = jy*alpha01*jM/M1+(1-jy)*alpha02*jM/M2;
	scalar gamma = 1+1.0/alpha;
	
	U(i,e*N_F+0) = rho;
	U(i,e*N_F+1) = rho*u;
	U(i,e*N_F+2) = rho*v;
	U(i,e*N_F+3) = p/(gamma-1)+ 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
	U(i,e*N_F+4) = rho/(gamma-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma-1);
#endif
      }
    }
  }
#endif
}

void init_dg_rarecon_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){


  // Initialize a rarefaction moving towards the left (or downwards)
  scalar A0 = 0.00183;
  scalar Lx = 0.089*2.0/3.0;
  scalar K = 6;
  scalar yinterface =-(K+1)*Lx;//-0.07; // first interface location
  scalar delta=0.005;   // The diffusion layer thickness
  scalar u0 = 0;
  scalar v0 = 0;
  scalar p0 = 1e5;

  // Top material (rarefaction initialized here)
  scalar rho01   = 5.494;//1.351;//5.494;//
  scalar gamma01 = 1.276;//;1.093;//1.276;//;
  scalar alpha01 = 1/(gamma01-1);
  scalar M01     = 146.05;//34.76;//146.05;//
  scalar c01     = sqrt(gamma01*p0/rho01);
  
  // Bottom material (material 2)
  scalar rho02   = 1.351;//5.494;//1.351;//
  scalar gamma02 = 1.093;//1.276;//
  scalar alpha02 = 1/(gamma02-1);
  scalar M02     = 34.76;//146.05;//34.76;// // molecular weight
  scalar c02     = sqrt(gamma02*p0/rho02);
    
  // Parameters to choose
  scalar strength = 0.73; // ratio pc/p0
  scalar L = K*Lx;
  printf("L = h/lambda=%f, strength=%f\n",K,strength);
  
  // contact velocity (velocity behind the rarefaction, positive)
  scalar vc = 2.0*c01/(gamma01-1) * (1 - pow(strength, (gamma01-1)/(2.0*gamma01)));   

  // time at with the rarefaction is of length L
  scalar ti = 2.0/(gamma01+1) * 1.0/vc * L;
  printf("interaction time with interface, ti = %f\n",ti);
  
  // solve for the origin of rarefaction
  scalar yR = yinterface+c01*ti;

  // Head and tail positions at t0
  scalar t0 = 0.80*ti; // alternatively: t0 = fudge * (yR-A0-yinterface)/c01; //0<fudge<1
  scalar yH0 = -c01*t0 + yR;    
  scalar yT0 = -(c01 - (gamma01+1)/2.0 * vc) * t0 + yR;

  // reflection coefficient
  scalar R = (1- rho01*c01/(rho02*c02))/(1+rho01*c01/(rho02*c02));
  scalar vcoord = -(1-R)*vc; // coordinate shift upwards
  
  scalar xc=0, yc=0, x=0, y=0, u, v, rho, p;
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){

#ifdef ONED
      A0 = 0;
      yc = XYZCen(e,0);
      y  = XYZNodes(i,e*D+0);
#elif TWOD
      xc = XYZCen(e,0);
      x  = XYZNodes(i,e*D+0);
      yc = XYZCen(e,1);
      y  = XYZNodes(i,e*D+1);
#endif

      // vertical distance from interface
      scalar d = ((delta+A0*sin(2*M_PI*x/Lx-M_PI/2))-y+yinterface)/(2*delta);
      
      // Calculate volume fractions
      scalar vol=0;
      if      ((d<1)&&(d>0)) vol = exp(log(1e-16)*pow(fabs(d),8));
      else if (d<=0)         vol = 1;
      else                   vol = 0;
      
      scalar jx  = 1-vol;
      scalar rho = jx*rho02+(1-jx)*rho01;
      scalar jy  = jx*rho02/(jx*rho02+(1-jx)*rho01);      // mass fraction
      scalar jM  = 1/(jy/M02+(1-jy)/M01);                 // total molecular weight
      
      scalar alpha = jy*alpha02*jM/M02+(1-jy)*alpha01*jM/M01;
      scalar gamma = 1+1.0/alpha;
      scalar c = sqrt(gamma*p0/rho);

      // Velocity, density and pressure modifications for rarefaction
      if (y<yH0){ // region in front of the rarefaction
	u = 0;
	v = 0;
	rho = rho;
	p = p0;	
      }
      else if (y>yT0){ // region behind the rarefaction
	u = 0;
	v = vc;
	rho = rho * pow(strength,1.0/gamma);//rho * pow(1 - (gamma-1)/2.0 * fabs(vc/c), 2.0/(gamma-1));
	p   = p0  * strength;
      }
      else{ // inside the rarefaction
	u = 0;
	v = 2.0/(gamma+1)*(c + (y-yR)/t0);
	rho = rho * pow(1 - (gamma-1)/2.0 * fabs(v/c), 2.0/(gamma-1));
	p   = p0  * pow(1 - (gamma-1)/2.0 * fabs(v/c), 2.0*gamma/(gamma-1));
      }
      v = v + vcoord;
      
#ifdef ONED
      U(i,e*N_F+0) = rho;
      U(i,e*N_F+1) = rho*v;
      U(i,e*N_F+2) = p/(gamma-1)+ 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
      U(i,e*N_F+3) = rho/(gamma-1);
#elif GAMNCON
      U(i,e*N_F+3) = 1.0/(gamma-1);
#endif
      // Mass fractions
      U(i,e*N_F+4) = (1-jx)*rho;
      
#elif TWOD
      U(i,e*N_F+0) = rho;
      U(i,e*N_F+1) = rho*u;
      U(i,e*N_F+2) = rho*v;
      U(i,e*N_F+3) = p/(gamma-1)+ 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
      U(i,e*N_F+4) = rho/(gamma-1);
#elif GAMNCON
      U(i,e*N_F+4) = 1.0/(gamma-1);
#endif
      // Mass fractions
      U(i,e*N_F+5) = (1-jx)*rho;
#endif
    }
  }
}



void init_dg_blastrm_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

  if ((N_F!=5)&&(N_F!=4)) printf("You are setting up the wrong problem. N_F =%i != 5 or 4.\n",N_F);

  // Initialize
  scalar A0 = 0.00183;
  scalar Lx = 0.089*2.0/3.0;
  scalar vcoord = 0; //-115.26;//72; // coordinate shift upwards
  scalar K = 6;
  scalar h = K*Lx;
  scalar yinterface =-h;//-0.07; // first interface location
  scalar delta=0.005;   // The diffusion layer thickness
  scalar u = 0;
  scalar v = 0+vcoord;
  scalar p = 1e5;

  // Top material: blast initialized here
  scalar rho01   = 1.351;//5.494;//
  scalar u01     = u;
  scalar v01     = v;  
  scalar gamma01 = 1.093;//1.276;//;
  scalar alpha01 = 1/(gamma01-1);
  scalar M01     = 34.76;//146.05;//

  // Bottom material (material 2)
  scalar rho02   = 5.494;//1.351;//
  scalar u02     = u;
  scalar v02     = v;  
  scalar gamma02 = 1.276;//1.093;//
  scalar alpha02 = 1/(gamma02-1);
  scalar M02     = 146.05;//34.76;//  // molecular weight
  
  // Explosion energy parameters
  scalar blstpos = 0.0;
  // scalar alpha = 2./3.; // = 2/3 for planar, 1/2 for cyl, 2/5 for sph.
  // scalar Q = 0.8119080; // for gamma = 1.4 //0.66927; // for alpha = 2/3 and gamma = 5/3
  // scalar Ms = 10;    // Mach number when blast front gets to interface
  // scalar ps = p*((2.0*gamma01*Ms*Ms-(gamma01-1))/(gamma01+1));  // pressure at shock in Pa
  // scalar xi = fabs(blstpos-yinterface); // distance btw blast initial position and interface
  // scalar Ex  = pow(Q,2.0/alpha)*0.5*(gamma01+1)*ps/(alpha*alpha)*pow(xi,2.0/alpha-2.0);

  scalar pratio = 500; // pressure ratio at blast = px/p0
  scalar Px = pratio*p;
  scalar Ex = Px/(gamma01-1);
  //scalar Dxx = 0.001; // energy initially deposited in Dxx
  
  printf("h/lambda=%f, pratio=%f, px=%e and Ex=%e\n",K,pratio,Px, Ex);
  
  scalar xc=0, yc=0, x=0, y=0;
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){

#ifdef ONED
      A0 = 0;
      yc = XYZCen(e,0);
      y  = XYZNodes(i,e*D+0);
#elif TWOD
      xc = XYZCen(e,0);
      x  = XYZNodes(i,e*D+0);
      yc = XYZCen(e,1);
      y  = XYZNodes(i,e*D+1);
#endif
      
      if(yc >= blstpos){ // blast region

#ifdef ONED
	U(i,e*N_F+0) = rho01;
	U(i,e*N_F+1) = rho01*v01;
	U(i,e*N_F+2) = Ex;///Dxx;
#ifdef GAMCONS
	U(i,e*N_F+3) = rho01/(gamma01-1);
#elif GAMNCON
	U(i,e*N_F+3) = 1.0/(gamma01-1);
#endif
	// Mass fractions
	U(i,e*N_F+4) = 1*rho01;//0;//1*rho01;
#elif TWOD
	U(i,e*N_F+0) = rho01;
	U(i,e*N_F+1) = rho01*u01;
	U(i,e*N_F+2) = rho01*v01;
	U(i,e*N_F+3) = Ex;///Dxx;
#ifdef GAMCONS
	U(i,e*N_F+4) = rho01/(gamma01-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma01-1);
#endif
	// Mass fractions
	U(i,e*N_F+5) = 1*rho01;//0;//1*rho01;
#endif
      }
      else{
	// vertical distance from interface
	scalar d = ((delta+A0*sin(2*M_PI*x/Lx-M_PI/2))-y+yinterface)/(2*delta);
	
	// Calculate volume fractions
	scalar vol=0;
	if      ((d<1)&&(d>0)) vol = exp(log(1e-16)*pow(fabs(d),8));
	else if (d<=0)         vol = 1;
	else                   vol = 0;

	scalar jx  = 1-vol;
	scalar rho = jx*rho02+(1-jx)*rho01;
	scalar jy  = jx*rho02/(jx*rho02+(1-jx)*rho01);      // mass fraction
	scalar jM  = 1/(jy/M02+(1-jy)/M01);                 // total molecular weight

	scalar alpha = jy*alpha02*jM/M02+(1-jy)*alpha01*jM/M01;
	scalar gamma = 1+1.0/alpha;
  
#ifdef ONED
	U(i,e*N_F+0) = rho;
	U(i,e*N_F+1) = rho*v;
	U(i,e*N_F+2) = p/(gamma-1)+ 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
	U(i,e*N_F+3) = rho/(gamma-1);
#elif GAMNCON
	U(i,e*N_F+3) = 1.0/(gamma-1);
#endif
	// Mass fractions
	U(i,e*N_F+4) = (1-jx)*rho;//(1-jx)*rho; // jy is mass fraction of 2
#elif TWOD
	U(i,e*N_F+0) = rho;
	U(i,e*N_F+1) = rho*u;
	U(i,e*N_F+2) = rho*v;
	U(i,e*N_F+3) = p/(gamma-1)+ 0.5*rho*(u*u+v*v);
#ifdef GAMCONS
	U(i,e*N_F+4) = rho/(gamma-1);
#elif GAMNCON
	U(i,e*N_F+4) = 1.0/(gamma-1);
#endif
	// Mass fractions
	U(i,e*N_F+5) = (1-jx)*rho;//(1-jx)*rho;
#endif
      }
    }
  }
}


void init_dg_sinephi_passive(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  scalar rho     = 1.0;
  scalar u       = 1.0;
  scalar gamma = constants::GLOBAL_GAMMA;
  scalar phi     = 0.0;
  scalar sinephix= 0.0; // sine perturbation on phi
  scalar sinephiy= 0.0; // sine perturbation on phi
  scalar sinerhox= 0.0; // sine perturbation on rho
  scalar sinerhoy= 0.0; // sine perturbation on rho
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
      scalar y = XYZNodes(i,e*D+1);
      for(int b = 0; b < N_S; b++) if ((xS(b,0) <= x) && (x < xS(b+1,0)))  ind = b;
      if ((x<-1)||(x>1)) {sinephix = 0; sinerhox=0;}
      else               {sinephix = Aphi*sin(M_PI*x); sinephiy = Aphi*sin(M_PI*y); sinerhox = Arho*sin(4*M_PI*x); sinerhoy = Arho*sin(4*M_PI*y);}
#ifdef ONED
      if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);
      S(0,0) = (rho+sinerhox);
      S(1,0) = (rho+sinerhox)*u;
      S(2,0) = 1.0/(gamma-1.0)*p + 0.5*(rho+sinerhox)*u*u;
      S(3,0) = (rho+sinerhox)*(phi+sinephix);
      S(4,0) = phi+sinephix;
#elif TWOD
      if (N_F!=6) printf("You are setting up the wrong problem. N_F =%i != 6.\n",N_F);
      scalar v = 1;
      u = 1;
      S(0,0) = (rho+sinerhox+sinerhoy);
      S(1,0) = (rho+sinerhox+sinerhoy)*u;
      S(2,0) = (rho+sinerhox+sinerhoy)*v;
      S(3,0) = 1.0/(gamma-1.0)*p + 0.5*(rho+sinerhox+sinerhoy)*(u*u+v*v);
      S(4,0) = (rho+sinerhox+sinerhoy)*(phi+sinephix+sinephiy);
      S(5,0) = phi+sinephix+sinephiy;
#endif
      for(int k = 0; k < N_F; k++) U(i,e*N_F+k) = S(k,ind);
    }
  }
}


void init_dg_sodmono_passive(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  scalar gamma = constants::GLOBAL_GAMMA;
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


void init_dg_stffrnt_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

  // Left state
  scalar rhoL  = 1.0;
  scalar uL    = 1.0;
  scalar gammaL= 1.4;
  scalar pinfL = 1.6;
  scalar pL    = 1.0;
  scalar EtL   = 1.0/(gammaL-1.0)*pL + gammaL*pinfL/(gammaL-1) + 0.5*rhoL*uL*uL;
  scalar GL    = 1.0/(gammaL-1.0);

  
  // Right state
  scalar rhoR   = 0.125;
  scalar uR     = 1.0;
  scalar gammaR = 1.6;
  scalar pinfR  = 0;
  scalar pR     = 1.0;
  scalar EtR    = 1.0/(gammaR-1.0)*pR + gammaR*pinfR/(gammaR-1)  + 0.5*rhoR*uR*uR;
  scalar GR     = 1.0/(gammaR-1.0);

  
  scalar xc=0;
  for(int e = 0; e < N_E; e++){
    xc = XYZCen(e,0);
    for(int i = 0; i < N_s; i++){
      if ((-0.5<xc)&&(xc<0.5)){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = EtL;
	U(i,e*N_F+3) = GL;
	U(i,e*N_F+4) = gammaL*pinfL/(gammaL-1);
      }
      else {
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = EtR;
	U(i,e*N_F+3) = GR;
	U(i,e*N_F+4) = gammaR*pinfR/(gammaR-1);
      }
    }
  }
}

void init_dg_stfshck_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

  // Initial condition taken from Cocchi 1996, Shyue JCP 1998
  // Also found in Eric's thesis p. 62
  
  // Left state
  scalar rhoL  = 1.241;
  scalar uL    = 0.0;
  scalar gammaL= 1.4;
  scalar pinfL = 0;
  scalar pL    = 2.753;
  scalar EtL   = 1.0/(gammaL-1.0)*pL + gammaL*pinfL/(gammaL-1) + 0.5*rhoL*uL*uL;
  scalar GL    = 1.0/(gammaL-1.0);

  
  // Right state
  scalar rhoR   = 0.991;
  scalar uR     = 0.0;
  scalar gammaR = 5.5;
  scalar pinfR  = 1.505;
  scalar pR     = 3.059*1e-4;
  scalar EtR    = 1.0/(gammaR-1.0)*pR + gammaR*pinfR/(gammaR-1)  + 0.5*rhoR*uR*uR;
  scalar GR     = 1.0/(gammaR-1.0);

  
  scalar xc=0;
  for(int e = 0; e < N_E; e++){
    xc = XYZCen(e,0);
    for(int i = 0; i < N_s; i++){
      if (xc<0.0){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = EtL;
	U(i,e*N_F+3) = GL;
	U(i,e*N_F+4) = gammaL*pinfL/(gammaL-1);
      }
      else {
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = EtR;
	U(i,e*N_F+3) = GR;
	U(i,e*N_F+4) = gammaR*pinfR/(gammaR-1);
      }
    }
  }
}

void init_dg_stfbubl_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U){

  // Shock-bubble interaction: air bubble in water
  // From: http://dx.doi.org/10.3850/978-981-07-2826-7_213
  // Pressure and pinf normalized by (rho_air*cs_air^2)

  // Material properties
  // air at 300K/27C from http://www.mhtl.uwaterloo.ca/old/onlinetools/airprop/airprop.html
  scalar rho_air = 1.1765;
  scalar gamma_air = 1.4;
  scalar patm = 101325;
  scalar cs_air = sqrt(gamma_air*patm/rho_air);

  // water at 300K
  scalar rho_water = 996; // use 970 if gelatine/water mixture as in Bourne1992
  scalar gamma_water = 5.5;// Shahab uses different EoS model and so 2.35;
  scalar pinf_water = 492115000;// Shahab uses: 1e9;
  scalar cs_water = sqrt(gamma_water*(patm+pinf_water)/rho_water);
  
  // Shock properties
  scalar ratio = 2000;//10000; // pressure ratio at shock
  printf("Pressure ratio = %g\n",ratio);
  
  // post-shock state
  scalar pSwater = ratio*patm;
  scalar rhoSwater = rho_water * ((gamma_water+1)/(gamma_water-1) * (pSwater+pinf_water)/(patm+pinf_water) + 1) / ((gamma_water+1)/(gamma_water-1) + (pSwater+pinf_water)/(patm+pinf_water));
  scalar uSwater = cs_water/gamma_water * (ratio-1)*patm/(patm+pinf_water) / sqrt((gamma_water+1)/(2*gamma_water) * (ratio-1) * patm/(patm+pinf_water) + 1);
  scalar rhoS  = rhoSwater/rho_air;
  scalar uS    = uSwater/cs_air;
  scalar vS    = 0;
  scalar gammaS= gamma_water;
  scalar pinfS = pinf_water/(rho_air*cs_air*cs_air);
  scalar pS    = pSwater/(rho_air*cs_air*cs_air);
  scalar EtS   = 1.0/(gammaS-1.0)*pS + gammaS*pinfS/(gammaS-1) + 0.5*rhoS*uS*uS;
  scalar GS    = 1.0/(gammaS-1.0);
  scalar xshckpos = -2;

  // background water
  scalar rhoW   = rho_water/rho_air;
  scalar uW     = 0.0;
  scalar vW     = 0.0;
  scalar gammaW = gamma_water;
  scalar pinfW  = pinf_water/(rho_air*cs_air*cs_air);
  scalar pW     = patm/(rho_air*cs_air*cs_air);
  scalar EtW    = 1.0/(gammaW-1.0)*pW + gammaW*pinfW/(gammaW-1)  + 0.5*rhoW*uW*uW;
  scalar GW     = 1.0/(gammaW-1.0);

  // Bubble properties
  scalar rhoB   = rho_air/rho_air;
  scalar uB     = 0.0;
  scalar vB     = 0.0;
  scalar gammaB = gamma_air;
  scalar pinfB  = 0;
  scalar pB     = patm/(rho_air*cs_air*cs_air);
  scalar EtB    = 1.0/(gammaB-1.0)*pB + gammaB*pinfB/(gammaB-1)  + 0.5*rhoB*uB*uB;
  scalar GB     = 1.0/(gammaB-1.0);
  scalar radius = 1;
 
  scalar xc=0,yc=0;
  for(int e = 0; e < N_E; e++){
    xc = XYZCen(e,0);
    yc = XYZCen(e,1);
    for(int i = 0; i < N_s; i++){
      if (xc*xc+yc*yc<radius*radius){ // inside the bubble
	U(i,e*N_F+0) = rhoB;
	U(i,e*N_F+1) = rhoB*uB;
	U(i,e*N_F+2) = rhoB*vB;
	U(i,e*N_F+3) = EtB;
	U(i,e*N_F+4) = GB;
	U(i,e*N_F+5) = gammaB*pinfB/(gammaB-1);
      }
      else if (xc<xshckpos){
	U(i,e*N_F+0) = rhoS;
	U(i,e*N_F+1) = rhoS*uS;
	U(i,e*N_F+2) = rhoS*vS;
	U(i,e*N_F+3) = EtS;
	U(i,e*N_F+4) = GS;
	U(i,e*N_F+5) = gammaS*pinfS/(gammaS-1);
      }
      else {
	U(i,e*N_F+0) = rhoW;
	U(i,e*N_F+1) = rhoW*uW;
	U(i,e*N_F+2) = rhoW*vW;
	U(i,e*N_F+3) = EtW;
	U(i,e*N_F+4) = GW;
	U(i,e*N_F+5) = gammaW*pinfW/(gammaW-1);
      }
    }
  }
}


void init_dg_euler1D_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

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

void init_dg_euler2D_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

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

void init_dg_sodtube_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

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


void init_dg_explode_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

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




void init_dg_brio_wu_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

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

void init_dg_alfvenw_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

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

void init_dg_fastshk_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

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

void init_dg_slowshk_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

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

void init_dg_fastrar_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

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

void init_dg_slowrar_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U){

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


void init_dg_ovortex_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

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

void init_dg_mhdroto_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U){

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

