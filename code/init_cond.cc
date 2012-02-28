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


void init_dg_simplew_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 8.\n",N_F);
  
  // Initial conditions (see CFD2 project 1)

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
	if      (model==0) U(i,e*N_F+3) = (scalar)rho*gamma;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gamma-1);
      }
      else if ((x>-1.5)&&(x<-0.5)){
	u = -1.0/gamma*(1-tanh((x+1)/(0.25-(x+1)*(x+1))));
	a = u*(gamma-1)/2.0 + 1.0;
	rho = gamma*pow(a,2.0/(gamma-1));;
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
	if      (model==0) U(i,e*N_F+3) = (scalar)rho*gamma;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gamma-1);
      }
      else if ((x>=-0.5)&&(x<=0.5)){
	u = 0;
	a = 0;
	rho = gamma;
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)1.0/(gamma-1);
	if      (model==0) U(i,e*N_F+3) = (scalar)rho*gamma;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gamma-1);
      }
      else if ((x>0.5)&&(x<1.5)){
	u = 1.0/gamma*(1+tanh((x-1)/(0.25-(x-1)*(x-1))));;
	a = 1 - (gamma-1)/2.0*u;
	rho = gamma*pow(a,2.0/(gamma-1));
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
	if      (model==0) U(i,e*N_F+3) = (scalar)rho*gamma;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gamma-1);
      }
      else if (x>=1.5){
	u = 2.0/gamma;
	a = 1 - (gamma-1)/gamma;
	rho = gamma*pow(a,2.0/(gamma-1));
	U(i,e*N_F+0) = (scalar)rho;
	U(i,e*N_F+1) = (scalar)rho*u;
	U(i,e*N_F+2) = (scalar)rho*a*a/(gamma*(gamma-1)) + 0.5*rho*u*rho*u/rho;
	if      (model==0) U(i,e*N_F+3) = (scalar)rho*gamma;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gamma-1);
      }
    }
  }
}

void init_dg_sodtube_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);
  
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
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if (x<0){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = EtL ;
	if      (model==0) U(i,e*N_F+3) = (scalar)rhoL*gammaL;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gammaL-1);
      }
      else if (x>=0){
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = EtR ;
	if      (model==0) U(i,e*N_F+3) = (scalar)rhoR*gammaR;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gammaR-1);
      }
    }
  }
}

void init_dg_contact_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);
  
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
  
  for(int e = 0; e < N_E; e++){
    scalar x = XYZNodes(0,e*D+0);
    for(int i = 0; i < N_s; i++){
      if (x<0){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = EtL ;
	if      (model==0) U(i,e*N_F+3) = (scalar)rhoL*gammaL;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gammaL-1);
      }
      else if (x>=0){
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = EtR ;
	if      (model==0) U(i,e*N_F+3) = (scalar)rhoR*gammaR;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gammaR-1);
      }
    }
  }
}


void init_dg_matfrnt_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);
  
  // Left state
  scalar rhoL  = 1.0;
  scalar uL    = 1.0;
  scalar gammaL= 1.2;
  scalar pL    = 1.0;
  scalar EtL   = 1.0/(gammaL-1.0)*pL + 0.5*rhoL*uL*uL;
    
  // Right state
  scalar rhoR   = 1.0;
  scalar uR     = 1.0;
  scalar gammaR = 1.6;
  scalar pR     = 1.0;
  scalar EtR    = 1.0/(gammaR-1.0)*pR + 0.5*rhoR*uR*uR;
  
  for(int e = 0; e < N_E; e++){
    scalar x = XYZNodes(0,e*D+0);
    for(int i = 0; i < N_s; i++){
      if (x<1E-6){
	U(i,e*N_F+0) = rhoL;
	U(i,e*N_F+1) = rhoL*uL;
	U(i,e*N_F+2) = EtL ;
	if      (model==0) U(i,e*N_F+3) = (scalar)rhoL*gammaL;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gammaL-1);

      }
      else if (x>=1E-6){
	U(i,e*N_F+0) = rhoR;
	U(i,e*N_F+1) = rhoR*uR;
	U(i,e*N_F+2) = EtR ;
	if      (model==0) U(i,e*N_F+3) = (scalar)rhoR*gammaR;
	else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gammaR-1);

      }
    }
  }
}

void init_dg_sinegam_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);

  scalar rho     = 1.0;
  scalar u       = 1.0;
  scalar gamma0  = 1.4;
  scalar sinegam = 0.0;  // sine perturbation on gamma
  scalar A       = 0.10; // amplitude of the perturbation
  scalar p       = 1.0;
  scalar Et      = 1.0/(gamma0+sinegam-1.0)*p + 0.5*rho*u*u;
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if ((x<-1)||(x>1)) sinegam = 0;
      else               sinegam = A*sin(M_PI*x);
      U(i,e*N_F+0) = rho;
      U(i,e*N_F+1) = rho*u;
      U(i,e*N_F+2) = 1.0/(gamma0+sinegam-1.0)*p + 0.5*rho*u*u;
      if      (model==0) U(i,e*N_F+3) = (scalar)rho*(gamma0+sinegam);
      else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gamma0+sinegam-1);
    }
  }
}
             
void init_dg_expogam_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  if (N_F!=4) printf("You are setting up the wrong problem. N_F =%i != 4.\n",N_F);

  scalar rho     = 1.0;
  scalar u       = 1.0;
  scalar gamma0  = 1.4;
  scalar expogam = 0.0;  // exp perturbation on gamma
  scalar A       = 0.10; // amplitude of the perturbation
  scalar p       = 1.0;
  scalar Et      = 1.0/(gamma0+expogam-1.0)*p + 0.5*rho*u*u;
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      expogam = 0.0*exp(-x*x/(0.2*0.2));
      U(i,e*N_F+0) = rho;
      U(i,e*N_F+1) = rho*u;
      U(i,e*N_F+2) = 1.0/(gamma0+expogam-1.0)*p + 0.5*rho*u*u;
      if      (model==0) U(i,e*N_F+3) = (scalar)rho*(gamma0+expogam);
      else if (model==1) U(i,e*N_F+3) = (scalar)1.0/(gamma0+expogam-1);
    }
  }
}

void init_dg_sinephi_passive(const int N_s, const int N_E, const int N_F, const int D, scalar &gamma, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U){

  if (N_F!=5) printf("You are setting up the wrong problem. N_F =%i != 5.\n",N_F);

  scalar rho     = 1.0;
  scalar u       = 1.0;
  gamma          = 1.4;
  scalar phi     = 0.0;
  scalar sinephi = 0.0;  // sine perturbation on phi
  scalar sinerho = 0.0;  // sine perturbation on rho
  scalar Aphi    = 1.0; // amplitude of the perturbation
  scalar Arho    = 0.10; // amplitude of the perturbation
  scalar p       = 1.0;
  scalar Et      = 1.0/(gamma-1.0)*p + 0.5*rho*u*u;
  
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      scalar x = XYZNodes(i,e*D+0);
      if ((x<-1)||(x>1)) {sinephi = 0; sinerho=0;}
      else               {sinephi = Aphi*sin(M_PI*x); sinerho = Arho*sin(4*M_PI*x);}
      U(i,e*N_F+0) = (rho+sinerho);
      U(i,e*N_F+1) = (rho+sinerho)*u;
      U(i,e*N_F+2) = 1.0/(gamma-1.0)*p + 0.5*(rho+sinerho)*u*u;
      U(i,e*N_F+3) = (rho+sinerho)*(phi+sinephi);
      U(i,e*N_F+4) = phi+sinephi;
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
