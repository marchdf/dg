#include <cpu_kernels.h>
#include <cstdlib>
#include <stdio.h>
#include <algorithm>

// Kernel definitions
void cpu_equal(int N_s, int N_E, int N_F, scalar* A, scalar* B){
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      for(int fc = 0; fc < N_F; fc++){
	A[(e*N_F+fc)*N_s+i] = B[(e*N_F+fc)*N_s+i];
      }
    }
  }
}

void cpu_add(int N_s, int N_E, int N_F, scalar* A, scalar* B, scalar c){

  // A = A + c*B
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      for(int fc = 0; fc < N_F; fc++){
	A[(e*N_F+fc)*N_s+i] =  A[(e*N_F+fc)*N_s+i] + c*B[(e*N_F+fc)*N_s+i];
      }
    }
  }
}


void cpu_mapToFace_shallow(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){

  for(int t = 0; t < M_T; t++){
    for(int j = 0; j < M_s; j++){
      for(int fc = 0; fc < N_F; fc++){
	for(int d = 0; d < 2; d++){
	  int idx = -1;
	  int face;
	  face= ((t*N_F+fc)*2+d)*M_s+j;
	  idx = map[face];
	  if(idx != -1){
	    UF[face] = U[idx];
	  }
	  else if (idx == -1){
	    if      (fc == 0) UF[((t*N_F+fc)*2+1)*M_s+j] = UF[((t*N_F+fc)*2+0)*M_s+j]; // eta
	    else if (fc == 1) UF[((t*N_F+fc)*2+1)*M_s+j] = -UF[((t*N_F+fc)*2+0)*M_s+j]; // ux
	    else if (fc == 2) UF[((t*N_F+fc)*2+1)*M_s+j] = -UF[((t*N_F+fc)*2+0)*M_s+j]; // uy
	  }
	}
      }
    }
  }
}

void cpu_mapToFace_mhd(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){

  for(int t = 0; t < M_T; t++){
    for(int j = 0; j < M_s; j++){
      for(int fc = 0; fc < N_F; fc++){
	for(int d = 0; d < 2; d++){
	  int idx = -1;
	  int face;
	  face= ((t*N_F+fc)*2+d)*M_s+j;
	  idx = map[face];
	  if(idx != -1){
	    UF[face] = U[idx];
	  }
	}
      }
    }
  }
}


void cpu_mapToFace_multifluid(int M_s, int M_T, int N_F, int N_s, int boundaryMap, scalar* U, scalar* UF){

  // All other boundaries
  for(int t = 1; t < M_T-1; t++){ 
    for(int fc = 0; fc < N_F; fc++){    
      UF[(t*N_F+fc)*2+0] = U[((t-1)*N_F+fc)*N_s+1];
      UF[(t*N_F+fc)*2+1] = U[(t*N_F+fc)*N_s+0];
    }
  }

  // Start and end boundaries
  for(int fc = 0; fc < N_F; fc++){
    UF[(0*N_F+fc)*2+1]       = U[(0*N_F+fc)*N_s+0];     
    UF[((M_T-1)*N_F+fc)*2+0] = U[((M_T-2)*N_F+fc)*N_s+1];
    if      (boundaryMap == 0){      //farfield
      UF[(0*N_F+fc)*2+0]        = UF[(0*N_F+fc)*2+1];
      UF[((M_T-1)*N_F+fc)*2+1]  = UF[((M_T-1)*N_F+fc)*2+0];
    }
    else if (boundaryMap == M_T-1){  //periodic
      UF[(0*N_F+fc)*2+0]        = UF[((M_T-1)*N_F+fc)*2+0];
      UF[((M_T-1)*N_F+fc)*2+1]  = UF[(0*N_F+fc)*2+1];
    }
  }
}

void cpu_mapToFace_passive(int M_s, int M_T, int N_F, int N_s, int boundaryMap, scalar* U, scalar* UF){

  // All other boundaries
  for(int t = 1; t < M_T-1; t++){ 
    for(int fc = 0; fc < N_F; fc++){    
      UF[(t*N_F+fc)*2+0] = U[((t-1)*N_F+fc)*N_s+1];
      UF[(t*N_F+fc)*2+1] = U[(t*N_F+fc)*N_s+0];
    }
  }

  // Start and end boundaries
  for(int fc = 0; fc < N_F; fc++){
    UF[(0*N_F+fc)*2+1]       = U[(0*N_F+fc)*N_s+0];     
    UF[((M_T-1)*N_F+fc)*2+0] = U[((M_T-2)*N_F+fc)*N_s+1];
    if      (boundaryMap == 0){      //farfield
      UF[(0*N_F+fc)*2+0]        = UF[(0*N_F+fc)*2+1];
      UF[((M_T-1)*N_F+fc)*2+1]  = UF[((M_T-1)*N_F+fc)*2+0];
    }
    else if (boundaryMap == M_T-1){  //periodic
      UF[(0*N_F+fc)*2+0]        = UF[((M_T-1)*N_F+fc)*2+0];
      UF[((M_T-1)*N_F+fc)*2+1]  = UF[(0*N_F+fc)*2+1];
    }
  }
}


void cpu_boundary(int M_s, int N_F, int M_B, int* boundaryMap, scalar* UF){

  for(int t = 0; t < M_B; t++){
    int t1 = boundaryMap[t*2+0];
    int t2 = boundaryMap[t*2+1];
    for(int j = 0; j < M_s; j++){
      for(int fc = 0; fc < N_F; fc++){
	UF[((t1*N_F+fc)*2+1)*M_s+j] = UF[((t2*N_F+fc)*2+0)*M_s+j]; 
      }
    }
  }
}

void cpu_mapToElement(int N_s, int N_E, int N_F, scalar* Q, scalar* q){

  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      Q[(e*N_F+fc)*N_s+0] = q[(e*N_F+fc)*2+1];
      Q[(e*N_F+fc)*N_s+1] = q[((e+1)*N_F+fc)*2+0];
    }
  }
}

void cpu_collocationU(int D, int N_G, int N_s, int N_E, int N_F, scalar* Ug, scalar* dUg, scalar* phi, scalar* dphi, scalar* U){
  scalar sol = 0;

  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
      for(int fc = 0; fc < N_F; fc++){

	sol = 0.0;
	for(int i = 0; i < N_s; i++){
	  sol += phi[i*N_G+g] * U[(e*N_F+fc)*N_s+i];
	}
	Ug[(e*N_F+fc)*N_G+g] = sol;

	sol = 0.0;
	for(int a = 0; a < D; a++){
	  for(int i = 0; i < N_s; i++){
	    sol += dphi[(i*N_G+g)*D+a] * U[(e*N_F+fc)*N_s+i];
	  }
	  dUg[((e*N_F+fc)*N_G+g)*D+a] = sol;
	  sol = 0.0;
	}
      }
    }
  }
}

void cpu_collocationUF(int M_G, int M_s, int M_T, int N_F, scalar* UgF, scalar* psi, scalar* UF){

  scalar sol = 0;
  for(int t = 0; t < M_T; t++){
    for(int g = 0; g < M_G; g++){
      for(int fc = 0; fc < N_F; fc++){
	sol = 0;
	for(int d = 0; d < 2; d++){
	  for(int j = 0; j < M_s; j++){
	    sol += psi[j*M_G+g] * UF[((t*N_F+fc)*2+d)*M_s+j];
	  }
	  UgF[((t*N_F+fc)*2+d)*M_G+g] = sol;
	  sol = 0.0;
	}
      }
    }
  }
}

void cpu_evaluate_sf_shallow(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0){

  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
      scalar eta  =  Ug[(e*N_F+0)*N_G+g];
      
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;
      
      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = H0*Ug[(e*N_F+1)*N_G+g]; // u_x
      f[((e*N_F+1)*N_G+g)*D+0] = G0*eta; // eta 
      f[((e*N_F+2)*N_G+g)*D+0] = 0; 
      
      // Flux derive par rapport a y
      f[((e*N_F+0)*N_G+g)*D+1] = H0*Ug[(e*N_F+2)*N_G+g]; // u_y
      f[((e*N_F+1)*N_G+g)*D+1] = 0;
      f[((e*N_F+2)*N_G+g)*D+1] = G0*eta; // eta
    }
  }
}

scalar cpu_flux1_mhd(scalar rho, scalar u){return rho*u;}                                        // for f0X, f0X, f0X
scalar cpu_flux2_mhd(scalar rho, scalar u, scalar Bx, scalar pbar){return rho*u*u-Bx*Bx+pbar;}   // for f1X, f2Y, f3Z
scalar cpu_flux3_mhd(scalar rho, scalar u, scalar v, scalar Bx, scalar By){return rho*u*v-Bx*By;}// for f1Y, f1Z, f2X, f2Z, f3X, f3Y
scalar cpu_flux4_mhd(scalar u, scalar By, scalar v, scalar Bx){return u*By-v*Bx;}                // for f4Y, f4Z, f5X, f5Z, f6X, f6Y
scalar cpu_flux5_mhd(scalar EtplusPbar, scalar u, scalar vdotB, scalar Bx) {return EtplusPbar*u - vdotB*Bx;} // for f7X, f7Y, f7Z

void cpu_evaluate_sf_mhd(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma){

  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
      scalar rho = Ug[(e*N_F+0)*N_G+g];   
      scalar u   = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar v   = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
      scalar Bx  = Ug[(e*N_F+3)*N_G+g];
      scalar By  = Ug[(e*N_F+4)*N_G+g];
      scalar Et  = Ug[(e*N_F+5)*N_G+g];
      scalar vdotv  = u*u+v*v;       
      scalar BdotB  = Bx*Bx+By*By;
      scalar vdotB  = u*Bx+v*By;
      scalar divB   = 0;
      for(int alpha = 0; alpha < D; alpha++){
	divB += dUg[((e*N_F+3)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+0)*D+alpha] + dUg[((e*N_F+4)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+1)*D+alpha];
      }
      scalar p = (gamma-1)*(Et - 0.5*(rho*vdotv+BdotB));
      scalar Pbar = p+0.5*BdotB;
      scalar EtPbar = Et+Pbar;

      // Powel source
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;//-divB*Bx;
      s[(e*N_F+2)*N_G+g] = 0;//-divB*By;
      s[(e*N_F+3)*N_G+g] = 0;//-divB*u;
      s[(e*N_F+4)*N_G+g] = 0;//-divB*v;
      s[(e*N_F+5)*N_G+g] = 0;//-divB*vdotB;
            
      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = cpu_flux1_mhd(rho,u);                //rho*u; 
      f[((e*N_F+1)*N_G+g)*D+0] = cpu_flux2_mhd(rho,u,Bx,Pbar);        //rho*u*u-Bx*Bx+Pbar; 
      f[((e*N_F+2)*N_G+g)*D+0] = cpu_flux3_mhd(rho,u,v,Bx,By);        //rho*u*v-Bx*By; 
      f[((e*N_F+3)*N_G+g)*D+0] = 0;                                   //0;
      f[((e*N_F+4)*N_G+g)*D+0] = cpu_flux4_mhd(u,By,v,Bx);            //u*By-v*Bx;
      f[((e*N_F+5)*N_G+g)*D+0] = cpu_flux5_mhd(EtPbar,u,vdotB,Bx);    //EtplusPbar*u-vdotB*Bx;
      
      // Flux derive par rapport a y
      f[((e*N_F+0)*N_G+g)*D+1] = cpu_flux1_mhd(rho,v);                //rho*v;
      f[((e*N_F+1)*N_G+g)*D+1] = cpu_flux3_mhd(rho,v,u,By,Bx);        //rho*v*u-By*Bx;
      f[((e*N_F+2)*N_G+g)*D+1] = cpu_flux2_mhd(rho,v,By,Pbar);        //rho*v*v-By*By+Pbar;
      f[((e*N_F+3)*N_G+g)*D+1] = cpu_flux4_mhd(v,Bx,u,By);            //v*Bx-u*By;
      f[((e*N_F+4)*N_G+g)*D+1] = 0;                                   //0;
      f[((e*N_F+5)*N_G+g)*D+1] = cpu_flux5_mhd(EtPbar,v,vdotB,By);    //EtplusPbar*v-vdotB*By;
    }
  }
}

scalar cpu_flux1_multifluid(scalar rho, scalar u){return rho*u;}                  
scalar cpu_flux2_multifluid(scalar rho, scalar u, scalar p){return rho*u*u+p;} 
scalar cpu_flux3_multifluid(scalar EtplusP, scalar u) {return EtplusP*u;}
scalar cpu_flux4_multifluid(scalar rho, scalar u, scalar gamma) {return rho*u*gamma;}
scalar cpu_flux5_multifluid(scalar u, scalar gamma) {return u/(gamma-1);}

void cpu_evaluate_sf_multifluid(int D, int N_G, int N_E, int N_F, int model, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){

  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
      scalar rho = Ug[(e*N_F+0)*N_G+g];   
      scalar u   = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar Et  = Ug[(e*N_F+2)*N_G+g];
      scalar gamma=0;
      if      (model==0) gamma=Ug[(e*N_F+3)*N_G+g]/rho;
      else if (model==1) gamma=1+1.0/Ug[(e*N_F+3)*N_G+g];
      scalar p = (gamma-1)*(Et - 0.5*rho*u*u);
      scalar EtplusP = Et + p;
      scalar dudx = (dUg[(e*N_F+1)*N_G+g]-u*dUg[(e*N_F+0)*N_G+g])/rho*invJac[e*N_G+g];
      scalar dalphadx = dUg[(e*N_F+3)*N_G+g]*invJac[e*N_G+g];
      //dudx += dUg[((e*N_F+1)*N_G+g)*D+alpha];//*invJac[((e*N_G+g)*D+0)*D+alpha];
      //scalar dudx = dUg[e*N_G+g]*invJac[e*N_G+g];

      //if ((e==50)&&(g==1)) printf("dudx=%12.11f  dalphadx=%12.11f dudxi=%12.11f  invJac=%f\n", dudx, dalphadx, dUg[(e*N_F+1)*N_G+g], invJac[e*N_G+g]);
      //if ((e==100)&&(g==2)) printf("dudx=%12.11f  dudxi=%12.11f  invJac=%f\n", dudx, dUg[e*N_G+g], invJac[e*N_G+g]);
	    
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;
      if      (model==0) s[(e*N_F+3)*N_G+g] = 0;
      else if (model==1) s[(e*N_F+3)*N_G+g] = -u*dalphadx;

      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = cpu_flux1_multifluid(rho,u);       
      f[((e*N_F+1)*N_G+g)*D+0] = cpu_flux2_multifluid(rho,u,p);      
      f[((e*N_F+2)*N_G+g)*D+0] = cpu_flux3_multifluid(EtplusP,u);   
      if      (model==0) f[((e*N_F+3)*N_G+g)*D+0] = 0;//cpu_flux4_multifluid(rho,u,gamma);
      else if (model==1) f[((e*N_F+3)*N_G+g)*D+0] = 0;
    }
  }
}

scalar cpu_flux1_passive(scalar rho, scalar u){return rho*u;}                  
scalar cpu_flux2_passive(scalar rho, scalar u, scalar p){return rho*u*u+p;} 
scalar cpu_flux3_passive(scalar EtplusP, scalar u) {return EtplusP*u;}
scalar cpu_flux4_passive(scalar rho, scalar u, scalar phi) {return rho*u*phi;}

void cpu_evaluate_sf_passive(int D, int N_G, int N_E, int N_F, scalar gamma, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){

  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar Et    = Ug[(e*N_F+2)*N_G+g];
      scalar phic  = Ug[(e*N_F+3)*N_G+g]/rho;
      scalar phinc = Ug[(e*N_F+4)*N_G+g];
      scalar p     = (gamma-1)*(Et - 0.5*rho*u*u);
      scalar EtplusP = Et + p;
      scalar dphincdx = dUg[(e*N_F+4)*N_G+g]*invJac[e*N_G+g];

      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;
      s[(e*N_F+3)*N_G+g] = 0;
      s[(e*N_F+4)*N_G+g] = -u*dphincdx;

      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = cpu_flux1_passive(rho,u);       
      f[((e*N_F+1)*N_G+g)*D+0] = cpu_flux2_passive(rho,u,p);      
      f[((e*N_F+2)*N_G+g)*D+0] = cpu_flux3_passive(EtplusP,u);   
      f[((e*N_F+3)*N_G+g)*D+0] = cpu_flux4_passive(rho,u,phic);   
      f[((e*N_F+4)*N_G+g)*D+0] = 0;
    }
  }
}

void cpu_evaluate_q_shallow(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals){

  for(int t = 0; t < M_T; t++){
    for(int g = 0; g < M_G; g++){
      scalar nx = normals[t*2+0];
      scalar ny = normals[t*2+1];
      scalar etaL= UgF[((t*N_F+0)*2+0)*M_G+g];
      scalar etaR= UgF[((t*N_F+0)*2+1)*M_G+g];
      scalar uLn = UgF[((t*N_F+1)*2+0)*M_G+g] * nx + UgF[((t*N_F+2)*2+0)*M_G+g] * ny;
      scalar uRn = UgF[((t*N_F+1)*2+1)*M_G+g] * nx + UgF[((t*N_F+2)*2+1)*M_G+g] * ny;
      
      scalar h0 = H0;
      scalar g0 = G0;
      
      // first equation
      scalar qL = -0.5*h0*(uLn + uRn + sqrt(g0/h0)*(etaL-etaR)); // Left
      q[((t*N_F+0)*2+0)*M_G+g] = qL;
      q[((t*N_F+0)*2+1)*M_G+g] = -qL;
      // second
      qL = -0.5*g0*nx*(etaL+etaR+sqrt(h0/g0)*(uLn-uRn)); // Left
      q[((t*N_F+1)*2+0)*M_G+g] = qL;
      q[((t*N_F+1)*2+1)*M_G+g] = -qL;
      // third
      qL = -0.5*g0*ny*(etaL+etaR+sqrt(h0/g0)*(uLn-uRn)); // Left
      q[((t*N_F+2)*2+0)*M_G+g] = qL;
      q[((t*N_F+2)*2+1)*M_G+g] = -qL;
    }
  }
}



void cpu_evaluate_q_mhd(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar gamma, scalar* normals){
  
  for(int t = 0; t < M_T; t++){
    for(int g = 0; g < M_G; g++){
      scalar nx = normals[t*2+0];
      scalar ny = normals[t*2+1];
      scalar rhoL= UgF[((t*N_F+0)*2+0)*M_G+g];
      scalar rhoR= UgF[((t*N_F+0)*2+1)*M_G+g];
      scalar uL  = UgF[((t*N_F+1)*2+0)*M_G+g]/rhoL;
      scalar uR  = UgF[((t*N_F+1)*2+1)*M_G+g]/rhoR;
      scalar vL  = UgF[((t*N_F+2)*2+0)*M_G+g]/rhoL;
      scalar vR  = UgF[((t*N_F+2)*2+1)*M_G+g]/rhoR;
      scalar BxL = UgF[((t*N_F+3)*2+0)*M_G+g];
      scalar BxR = UgF[((t*N_F+3)*2+1)*M_G+g];
      scalar ByL = UgF[((t*N_F+4)*2+0)*M_G+g];
      scalar ByR = UgF[((t*N_F+4)*2+1)*M_G+g];
      scalar EtL = UgF[((t*N_F+5)*2+0)*M_G+g];
      scalar EtR = UgF[((t*N_F+5)*2+1)*M_G+g];
      scalar vdotvL = uL*uL+vL*vL;
      scalar vdotvR = uR*uR+vR*vR;
      scalar BdotBL = BxL*BxL+ByL*ByL;
      scalar BdotBR = BxR*BxR+ByR*ByR;
      scalar vdotBL = uL*BxL+vL*ByL;
      scalar vdotBR = uR*BxR+vR*ByR;
      scalar vdotnL = uL*nx+vL*ny;
      scalar vdotnR = uR*nx+vR*ny;
      scalar BdotnL = BxL*nx+ByL*ny;
      scalar BdotnR = BxR*nx+ByR*ny;
      scalar aveBdotn = 0.5*(BdotnL + BdotnR);      
      scalar pL = (gamma-1)*(EtL - 0.5*(rhoL*vdotvL+BdotBL));
      scalar pR = (gamma-1)*(EtR - 0.5*(rhoR*vdotvR+BdotBR));
      scalar PbarL = pL+0.5*BdotBL;
      scalar PbarR = pR+0.5*BdotBR;
      scalar EtPbarL = EtL+PbarL;
      scalar EtPbarR = EtR+PbarR;

      // Evaluate the right and left eigenvalues
      int sizevap = 8;
      scalar* vap = new scalar[2*sizevap];
      scalar alfenL  = BdotnL/sqrt(rhoL);
      scalar alfenR  = BdotnR/sqrt(rhoR);
      scalar a2L = (gamma*pL+BdotBL)/rhoL;
      scalar a2R = (gamma*pR+BdotBR)/rhoR;
      scalar cfL = sqrt(0.5*(a2L + sqrt( a2L*a2L - 4.0*gamma*pL*BdotnL*BdotnL/(rhoL*rhoL))));
      scalar cfR = sqrt(0.5*(a2R + sqrt( a2R*a2R - 4.0*gamma*pR*BdotnR*BdotnR/(rhoR*rhoR))));
      scalar csL = sqrt(0.5*(a2L - sqrt( a2L*a2L - 4.0*gamma*pL*BdotnL*BdotnL/(rhoL*rhoL))));
      scalar csR = sqrt(0.5*(a2R - sqrt( a2R*a2R - 4.0*gamma*pR*BdotnR*BdotnR/(rhoR*rhoR))));

      vap[0*sizevap+0] = fabs(vdotnL);
      vap[0*sizevap+1] = fabs(vdotnL) + alfenL;
      vap[0*sizevap+2] = fabs(vdotnL) - alfenL;
      vap[0*sizevap+3] = fabs(vdotnL) + cfL;
      vap[0*sizevap+4] = fabs(vdotnL) - cfL;
      vap[0*sizevap+5] = fabs(vdotnL) + csL;
      vap[0*sizevap+6] = fabs(vdotnL) - csL;
      vap[0*sizevap+7] = fabs(vdotnL);
      
      vap[1*sizevap+0] = fabs(vdotnR);
      vap[1*sizevap+1] = fabs(vdotnR) + alfenR;
      vap[1*sizevap+2] = fabs(vdotnR) - alfenR;
      vap[1*sizevap+3] = fabs(vdotnR) + cfR;
      vap[1*sizevap+4] = fabs(vdotnR) - cfR;
      vap[1*sizevap+5] = fabs(vdotnR) + csR;
      vap[1*sizevap+6] = fabs(vdotnR) - csR;
      vap[1*sizevap+7] = fabs(vdotnR);

      scalar maxvap = 0;
      for (int k = 0; k < 2*sizevap; k++){
	if (maxvap<vap[k]) maxvap = vap[k];
      }
      delete[] vap;

      // Upwinding on the source term
      // scalar upBdotnL = aveBdotn;
      // scalar upBdotnR = aveBdotn;
      scalar upBdotnL = 0.0;
      scalar upBdotnR = 0.0;
      if      (vdotnL >= 0) upBdotnL = BdotnL;
      else if (vdotnL <  0) upBdotnL = BdotnR;
      if      (vdotnR >= 0) upBdotnR = BdotnL;
      else if (vdotnR <  0) upBdotnR = BdotnR;  
      
      //
      // Evaluate the fluxes on the right and left
      //

      //first: fx = rho*u; fy = rho*v; fz = rho*w; 
      scalar qL = -0.5*((cpu_flux1_mhd(rhoL,uL) + cpu_flux1_mhd(rhoR,uR))*nx +
			(cpu_flux1_mhd(rhoL,vL) + cpu_flux1_mhd(rhoR,vR))*ny 
			-maxvap*(rhoR-rhoL));
      q[((t*N_F+0)*2+0)*M_G+g] = qL;
      q[((t*N_F+0)*2+1)*M_G+g] = -qL;
      
      //second: fx = rho*u*u+Bx*Bx+Pbar; fy = rho*v*u-By*Bx; fz = rho*w*u-Bz*Bx;
      qL = -0.5*((cpu_flux2_mhd(rhoL,uL,BxL,PbarL)  + cpu_flux2_mhd(rhoR,uR,BxR,PbarR) )*nx +
		 (cpu_flux3_mhd(rhoL,vL,uL,ByL,BxL) + cpu_flux3_mhd(rhoR,vR,uR,ByR,BxR))*ny
		 -maxvap*(rhoR*uR-rhoL*uL));
      q[((t*N_F+1)*2+0)*M_G+g] = qL  - BxL*upBdotnL;
      q[((t*N_F+1)*2+1)*M_G+g] = -qL + BxR*upBdotnR;

      //third: fx = rho*u*v-Bx*By; fy = rho*v*v-By*By+Pbar; fz = rho*w*v-Bz*By;
      qL = -0.5*((cpu_flux3_mhd(rhoL,uL,vL,BxL,ByL) + cpu_flux3_mhd(rhoR,uR,vR,BxR,ByR))*nx +
		 (cpu_flux2_mhd(rhoL,vL,ByL,PbarL)  + cpu_flux2_mhd(rhoR,vR,ByR,PbarR) )*ny
		 -maxvap*(rhoR*vR-rhoL*vL));
      q[((t*N_F+2)*2+0)*M_G+g] = qL  - ByL*upBdotnL;
      q[((t*N_F+2)*2+1)*M_G+g] = -qL + ByR*upBdotnR;
      
      //fourth: fx = 0;  fy = v*Bx-u*By; fz = w*Bx-u*Bz;
      qL = -0.5*((0                            + 0                           )*nx +
		 (cpu_flux4_mhd(vL,BxL,uL,ByL) + cpu_flux4_mhd(vR,BxR,uR,ByR))*ny
		 -maxvap*(BxR-BxL));
      q[((t*N_F+3)*2+0)*M_G+g] = qL  - uL*upBdotnL;
      q[((t*N_F+3)*2+1)*M_G+g] = -qL + uR*upBdotnR;

      //fifth: fx = u*By-v*Bx;  fy = 0; fz = w*By-v*Bz;
      qL = -0.5*((cpu_flux4_mhd(uL,ByL,vL,BxL) + cpu_flux4_mhd(uR,ByR,vR,BxR))*nx +
		 (0                            + 0                           )*ny 
		 -maxvap*(ByR-ByL));
      q[((t*N_F+4)*2+0)*M_G+g] = qL  - vL*upBdotnL;
      q[((t*N_F+4)*2+1)*M_G+g] = -qL + vR*upBdotnR;

      //sixth: fx = EtplusPbar*u-vdotB*Bx; fy = EtplusPbar*v-vdotB*By; fz = EtplusPbar*w-vdotB*Bz;
      qL = -0.5*((cpu_flux5_mhd(EtPbarL,uL,vdotBL,BxL) + cpu_flux5_mhd(EtPbarR,uR,vdotBR,BxR))*nx +
		 (cpu_flux5_mhd(EtPbarL,vL,vdotBL,ByL) + cpu_flux5_mhd(EtPbarR,vR,vdotBR,ByR))*ny 
		 -maxvap*(EtR-EtL));
      q[((t*N_F+5)*2+0)*M_G+g] = qL  - vdotBL*upBdotnL;
      q[((t*N_F+5)*2+1)*M_G+g] = -qL + vdotBR*upBdotnR; 
    }
  }
}

scalar cpu_fhll_multifluid(scalar UL, scalar SL, scalar FL, scalar UR, scalar SR, scalar FR){
  return (SR*FL-SL*FR+SL*SR*(UR-UL))/(SR-SL);
}
void cpu_evaluate_q_multifluid(int M_G, int M_T, int N_F, int flux, int model, scalar* q, scalar* UgF){
  
  for(int t = 0; t < M_T; t++){
    scalar rhoL= UgF[(t*N_F+0)*2+0];
    scalar rhoR= UgF[(t*N_F+0)*2+1];
    scalar uL  = UgF[(t*N_F+1)*2+0]/rhoL;
    scalar uR  = UgF[(t*N_F+1)*2+1]/rhoR;
    scalar EtL = UgF[(t*N_F+2)*2+0];
    scalar EtR = UgF[(t*N_F+2)*2+1];
    scalar gammaL = 0; scalar gammaR = 0;
    scalar alphaL = UgF[(t*N_F+3)*2+0]; scalar alphaR = UgF[(t*N_F+3)*2+1];
    if (model==0){
      gammaL = UgF[(t*N_F+3)*2+0]/rhoL;
      gammaR = UgF[(t*N_F+3)*2+1]/rhoR;}
    else if (model==1){
      gammaL = 1.0+1.0/alphaL;
      gammaR = 1.0+1.0/alphaR;}
    
    //printf("%f and %f\n",gammaL, gammaR);
    scalar pL = (gammaL-1)*(EtL - 0.5*rhoL*uL*uL);
    scalar pR = (gammaR-1)*(EtR - 0.5*rhoR*uR*uR);
    scalar EtPL = EtL+pL;
    scalar EtPR = EtR+pR;

    // Evaluate the right and left eigenvalues
    int sizevap = 4;
    scalar* vap = new scalar[2*sizevap];
    scalar aL = sqrt((gammaL*pL)/rhoL);
    scalar aR = sqrt((gammaR*pR)/rhoR);

    vap[0*sizevap+0] = fabs(uL) - aL;
    vap[0*sizevap+1] = fabs(uL);
    vap[0*sizevap+2] = fabs(uL);
    vap[0*sizevap+3] = fabs(uL) + aL;

    vap[1*sizevap+0] = fabs(uR) - aR;
    vap[1*sizevap+1] = fabs(uR);
    vap[1*sizevap+2] = fabs(uR);
    vap[1*sizevap+3] = fabs(uR) + aR;

    scalar maxvap = 0;
    for (int k = 0; k < 2*sizevap; k++){
      if (maxvap<vap[k]) maxvap = vap[k];
    }
    // scalar SL = std::min(vap[0*sizevap+0],vap[1*sizevap+0]);
    // scalar SR = std::max(vap[0*sizevap+3],vap[1*sizevap+3]);
    scalar SL = std::min(fabs(uL)-aL,fabs(uR)-aR);
    scalar SR = std::max(fabs(uL)+aR,fabs(uR)+aR);
    delete[] vap;
    
    //
    // Evaluate the fluxes on the right and left
    //

    // Local Lax-Freidrichs flux
    if (flux == 0){
      //first: fx = rho*u; 
      scalar qL = -0.5*(cpu_flux1_multifluid(rhoL,uL) + cpu_flux1_multifluid(rhoR,uR)
      			-maxvap*(rhoR-rhoL));
      q[(t*N_F+0)*2+0] = qL;
      q[(t*N_F+0)*2+1] = -qL;
      
      //second: fx = rho*u*u+Bx*Bx+Pbar; 
      qL = -0.5*(cpu_flux2_multifluid(rhoL,uL,pL)  + cpu_flux2_multifluid(rhoR,uR,pR)
      		 -maxvap*(rhoR*uR-rhoL*uL));
      q[(t*N_F+1)*2+0] = qL;
      q[(t*N_F+1)*2+1] = -qL;

      //third: fx = EtplusP*u; 
      qL = -0.5*(cpu_flux3_multifluid(EtPL,uL) + cpu_flux3_multifluid(EtPR,uR)
      		 -maxvap*(EtR-EtL));
      q[(t*N_F+2)*2+0] = qL; 
      q[(t*N_F+2)*2+1] = -qL;

      //fourth: 
      if      (model==0){ //fx = rho*u*gamma; 
	qL = -0.5*(cpu_flux4_multifluid(rhoL,uL,gammaL) + cpu_flux4_multifluid(rhoR,uR,gammaR)
		   -maxvap*(rhoR*gammaR-rhoL*gammaL));}
      else if (model==1){ 
	qL = -0.5*maxvap*(alphaL-alphaR);}
      q[(t*N_F+3)*2+0] = qL + 0.5*0.5*(uL+uR)*(alphaL-alphaR); 
      q[(t*N_F+3)*2+1] = -qL+ 0.5*0.5*(uL+uR)*(alphaL-alphaR);
      //strictly equiv:
      //q[(t*N_F+3)*2+0] = 0.5*(alphaL-alphaR)*(0.5*(uL+uR)-maxvap); 
      //q[(t*N_F+3)*2+1] = 0.5*(alphaL-alphaR)*(0.5*(uL+uR)+maxvap);
      
    }

    // Non-conservative flux
    else if (flux == 1){
      scalar pnc1=0, pnc2=0, pnc3=0, pnc4=0;
      scalar vnc4 = -0.5*(uL+uR)*(alphaL-alphaR);

      // define the flux
      if (SL > 0){
	pnc1 = cpu_flux1_multifluid(rhoL,uL);
	pnc2 = cpu_flux2_multifluid(rhoL,uL,pL);
	pnc3 = cpu_flux3_multifluid(EtPL,uL);
	pnc4 = -0.5*vnc4;
      }
      else if ((SL < 0)&&(SR > 0)){
	pnc1 = cpu_fhll_multifluid(rhoL, SL, cpu_flux1_multifluid(rhoL,uL),   rhoR, SR, cpu_flux1_multifluid(rhoR,uR));
	pnc2 = cpu_fhll_multifluid(  uL, SL, cpu_flux2_multifluid(rhoL,uL,pL),  uR, SR, cpu_flux2_multifluid(rhoR,uR,pR));
	pnc3 = cpu_fhll_multifluid( EtL, SL, cpu_flux3_multifluid(EtPL,uL),    EtR, SR, cpu_flux3_multifluid(EtPR,uR));
	pnc4 = cpu_fhll_multifluid(alphaL, SL, 0, alphaR, SR, 0) - 0.5*(SR+SL)/(SR-SL)*vnc4;
      }
      else if (SR < 0){
	pnc1 = cpu_flux1_multifluid(rhoR,uR);
	pnc2 = cpu_flux2_multifluid(rhoR,uR,pR);
	pnc3 = cpu_flux3_multifluid(EtPR,uR);
	pnc4 = 0.5*vnc4;
      }

      // first
      scalar qL = -pnc1;
      q[(t*N_F+0)*2+0] = qL;
      q[(t*N_F+0)*2+1] = -qL;
      
      //second: fx = rho*u*u+Bx*Bx+Pbar; 
      qL = -pnc2;
      q[(t*N_F+1)*2+0] = qL;
      q[(t*N_F+1)*2+1] = -qL;

      //third: fx = EtplusP*u; 
      qL = -pnc3;
      q[(t*N_F+2)*2+0] = qL; 
      q[(t*N_F+2)*2+1] = -qL;

      //fourth: 
      if      (model==0){ //fx = rho*u*gamma; 
	qL = -0.5*(cpu_flux4_multifluid(rhoL,uL,gammaL) + cpu_flux4_multifluid(rhoR,uR,gammaR)
		   -maxvap*(rhoR*gammaR-rhoL*gammaL));}
      else if (model==1){ //fx = u/(gamma-1);
      	qL = -pnc4;}
      q[(t*N_F+3)*2+0] = qL - 0.5*vnc4; 
      q[(t*N_F+3)*2+1] = -qL- 0.5*vnc4;
    }

    // Non-conservative Roe flux
    else if (flux == 2){
      scalar rhoRoe = sqrt(rhoL*rhoR);
      scalar uRoe = (sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar HL = (EtL + pL)/rhoL;
      scalar HR = (EtR + pR)/rhoR;
      scalar HRoe = (sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar alphaRoe = (sqrt(rhoL)*alphaL+sqrt(rhoR)*alphaR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar gammaRoe = 1+1.0/alphaRoe;
      scalar iL = pL*alphaL;
      scalar iR = pR*alphaR;
      scalar iRoe = (sqrt(rhoL)*iL+sqrt(rhoR)*iR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar aRoe = sqrt((gammaRoe-1)*(HRoe-0.5*uRoe*uRoe));
      scalar DpRoe= (gammaRoe-1)*(gammaRoe-1)*(alphaRoe*(iR-iL) - iRoe*(alphaR-alphaL));
      scalar pRoe = (gammaRoe-1)*iRoe;

      scalar* vapRoe = new scalar[4];
      vapRoe[0] = uRoe-aRoe;
      vapRoe[1] = uRoe;
      vapRoe[2] = uRoe+aRoe;
      vapRoe[3] = uRoe;

      scalar* vep = new scalar[4*4];
      vep[0*4+0] = 1;
      vep[0*4+1] = uRoe-aRoe;
      vep[0*4+2] = HRoe-uRoe*aRoe;
      vep[0*4+3] = 0;

      vep[1*4+0] = 1;
      vep[1*4+1] = uRoe;
      vep[1*4+2] = 0.5*uRoe*uRoe;
      vep[1*4+3] = 0;

      vep[2*4+0] = 1;
      vep[2*4+1] = uRoe+aRoe;
      vep[2*4+2] = HRoe+uRoe*aRoe;
      vep[2*4+3] = 0;
      
      vep[3*4+0] = 0;
      vep[3*4+1] = 0;
      vep[3*4+2] = pRoe;
      vep[3*4+3] = 1;
      
      scalar* aiRoe= new scalar[4];
      aiRoe[1] = (rhoR-rhoL) - DpRoe/(aRoe*aRoe);
      aiRoe[2] = (DpRoe + rhoRoe*aRoe*(uR-uL))/(2*aRoe*aRoe);
      aiRoe[0] = (DpRoe - rhoRoe*aRoe*(uR-uL))/(2*aRoe*aRoe);
      aiRoe[3] = alphaR-alphaL;
      // aiRoe[1]= (gammaRoe-1)/(aRoe*aRoe)*((rhoR-rhoL)*(HRoe-uRoe*uRoe)+uRoe*(rhoR*uR-rhoL*uL)-(EtR-EtL));
      // aiRoe[0]=1/(2*aRoe)*((rhoR-rhoL)*(uRoe+aRoe)-(rhoR*uR-rhoL*uL)-aRoe*aiRoe[1]);
      // aiRoe[2]=(rhoR-rhoL)-(aiRoe[0]+aiRoe[1]);

      //first: fx = rho*u; 
      scalar qL = 0;
      if      (uRoe>0)  qL = cpu_flux1_multifluid(rhoL,uL) + aiRoe[0]*vapRoe[0]*vep[0*4+0];
      else if (uRoe<=0) qL = cpu_flux1_multifluid(rhoR,uR) - aiRoe[2]*vapRoe[2]*vep[2*4+0];
      q[(t*N_F+0)*2+0] = -qL;
      q[(t*N_F+0)*2+1] = qL;

      //second: fx = rho*u*u+Bx*Bx+Pbar; 
      if      (uRoe>0)  qL = cpu_flux2_multifluid(rhoL,uL,pL) + aiRoe[0]*vapRoe[0]*vep[0*4+1];
      else if (uRoe<=0) qL = cpu_flux2_multifluid(rhoR,uR,pR) - aiRoe[2]*vapRoe[2]*vep[2*4+1];
      q[(t*N_F+1)*2+0] = -qL;
      q[(t*N_F+1)*2+1] = qL;

      //third: fx = EtplusP*u; 
      if      (uRoe>0)  qL = cpu_flux3_multifluid(EtPL,uL) + aiRoe[0]*vapRoe[0]*vep[0*4+2];
      else if (uRoe<=0) qL = cpu_flux3_multifluid(EtPR,uR) - aiRoe[2]*vapRoe[2]*vep[2*4+2];
      q[(t*N_F+2)*2+0] = -qL; 
      q[(t*N_F+2)*2+1] = qL;

      //fourth: 
      scalar alphastar = 0;
      scalar ustar = 0;
      scalar gstar = 0;
      if      (model==0){ //fx = rho*u*gamma; 
	qL = -0.5*(cpu_flux4_multifluid(rhoL,uL,gammaL) + cpu_flux4_multifluid(rhoR,uR,gammaR)
		   -maxvap*(rhoR*gammaR-rhoL*gammaL));}
      else if (model==1){ 
	if      (uRoe>0){
	  qL =  aiRoe[0]*vapRoe[0]*vep[0*4+3];
	  alphastar = alphaL + aiRoe[0]*vep[0*4+3];
	  ustar  = uL + aiRoe[0]*vep[0*4+1];
	  gstar = -0.5*(alphaR-alphaL)*(uR+uL);
	}
	else if (uRoe<=0){
	  qL = -aiRoe[2]*vapRoe[2]*vep[2*4+3];
	  alphastar = alphaR - aiRoe[2]*vep[2*4+3];
	  ustar  = uR + aiRoe[2]*vep[2*4+1];
	  gstar = 0.5*(alphaR-alphastar)*(uR+ustar);
	}
      }
      q[(t*N_F+3)*2+0] = 0;//qL + 0.5*gstar;
      q[(t*N_F+3)*2+1] = 0;//-qL  + 0.5*gstar;
     
      
      delete[] vapRoe;
      delete[] vep;
      delete[] aiRoe;
    }
  }
}

scalar cpu_fhll_passive(scalar UL, scalar SL, scalar FL, scalar UR, scalar SR, scalar FR){
  return (SR*FL-SL*FR+SL*SR*(UR-UL))/(SR-SL);
}
void cpu_evaluate_q_passive(int M_G, int M_T, int N_F, int flux, scalar gamma, scalar* q, scalar* UgF){
  
  for(int t = 0; t < M_T; t++){
    scalar rhoL   = UgF[(t*N_F+0)*2+0];
    scalar rhoR   = UgF[(t*N_F+0)*2+1];
    scalar uL     = UgF[(t*N_F+1)*2+0]/rhoL;
    scalar uR     = UgF[(t*N_F+1)*2+1]/rhoR;
    scalar EtL    = UgF[(t*N_F+2)*2+0];
    scalar EtR    = UgF[(t*N_F+2)*2+1];
    scalar phicL  = UgF[(t*N_F+3)*2+0]/rhoL;
    scalar phicR  = UgF[(t*N_F+3)*2+1]/rhoR;
    scalar phincL = UgF[(t*N_F+4)*2+0];
    scalar phincR = UgF[(t*N_F+4)*2+1];
    scalar pL     = (gamma-1)*(EtL - 0.5*rhoL*uL*uL);
    scalar pR     = (gamma-1)*(EtR - 0.5*rhoR*uR*uR);
    scalar EtPL   = EtL+pL;
    scalar EtPR   = EtR+pR;

    // Evaluate the right and left eigenvalues
    int sizevap = 4;
    scalar* vap = new scalar[2*sizevap];
    scalar aL = sqrt((gamma*pL)/rhoL);
    scalar aR = sqrt((gamma*pR)/rhoR);

    vap[0*sizevap+0] = fabs(uL) - aL;
    vap[0*sizevap+1] = fabs(uL);
    vap[0*sizevap+2] = fabs(uL);
    vap[0*sizevap+3] = fabs(uL) + aL;

    vap[1*sizevap+0] = fabs(uR) - aR;
    vap[1*sizevap+1] = fabs(uR);
    vap[1*sizevap+2] = fabs(uR);
    vap[1*sizevap+3] = fabs(uR) + aR;

    scalar maxvap = 0;
    for (int k = 0; k < 2*sizevap; k++){
      if (maxvap<vap[k]) maxvap = vap[k];
    }
    // scalar SL = std::min(vap[0*sizevap+0],vap[1*sizevap+0]);
    // scalar SR = std::max(vap[0*sizevap+3],vap[1*sizevap+3]);
    scalar SL = std::min(fabs(uL)-aL,fabs(uR)-aR);
    scalar SR = std::max(fabs(uL)+aR,fabs(uR)+aR);
    delete[] vap;
    
    //
    // Evaluate the fluxes on the right and left
    //

    // Local Lax-Freidrichs flux
    if (flux == 0){
      //first: fx = rho*u; 
      scalar qL = -0.5*(cpu_flux1_passive(rhoL,uL) + cpu_flux1_passive(rhoR,uR)
      			-maxvap*(rhoR-rhoL));
      q[(t*N_F+0)*2+0] = qL;
      q[(t*N_F+0)*2+1] = -qL;
      
      //second: fx = rho*u*u+Bx*Bx+Pbar; 
      qL = -0.5*(cpu_flux2_passive(rhoL,uL,pL)  + cpu_flux2_passive(rhoR,uR,pR)
      		 -maxvap*(rhoR*uR-rhoL*uL));
      q[(t*N_F+1)*2+0] = qL;
      q[(t*N_F+1)*2+1] = -qL;

      //third: fx = EtplusP*u; 
      qL = -0.5*(cpu_flux3_passive(EtPL,uL) + cpu_flux3_passive(EtPR,uR)
      		 -maxvap*(EtR-EtL));
      q[(t*N_F+2)*2+0] = qL; 
      q[(t*N_F+2)*2+1] = -qL;

      //fourth: fx = rho*u*phi;
      qL = -0.5*(cpu_flux4_passive(rhoL,uL,phicL) + cpu_flux4_passive(rhoR,uR,phicR)
      		 -maxvap*(phicR-phicL));
      q[(t*N_F+3)*2+0] = qL; 
      q[(t*N_F+3)*2+1] = -qL;

      //fifth:
      qL = -0.5*maxvap*(phincL-phincR);
      q[(t*N_F+4)*2+0] = qL + 0.5*0.5*(uL+uR)*(phincL-phincR); 
      q[(t*N_F+4)*2+1] = -qL+ 0.5*0.5*(uL+uR)*(phincL-phincR);
      //strictly equiv:
      //q[(t*N_F+3)*2+0] = 0.5*(alphaL-alphaR)*(0.5*(uL+uR)-maxvap); 
      //q[(t*N_F+3)*2+1] = 0.5*(alphaL-alphaR)*(0.5*(uL+uR)+maxvap);
    }

    // Non-conservative flux
    else if (flux == 1){
      scalar pnc1=0, pnc2=0, pnc3=0, pnc4=0, pnc5=0;
      scalar vnc5 = -0.5*(uL+uR)*(phincL-phincR);

      // define the flux
      if (SL > 0){
	pnc1 = cpu_flux1_passive(rhoL,uL);
	pnc2 = cpu_flux2_passive(rhoL,uL,pL);
	pnc3 = cpu_flux3_passive(EtPL,uL);
	pnc4 = cpu_flux4_passive(rhoL,uL,phicL);
	pnc5 = -0.5*vnc5;
      }
      else if ((SL < 0)&&(SR > 0)){
	pnc1 = cpu_fhll_passive( rhoL, SL, cpu_flux1_passive(rhoL,uL),        rhoR, SR, cpu_flux1_passive(rhoR,uR));
	pnc2 = cpu_fhll_passive(   uL, SL, cpu_flux2_passive(rhoL,uL,pL),       uR, SR, cpu_flux2_passive(rhoR,uR,pR));
	pnc3 = cpu_fhll_passive(  EtL, SL, cpu_flux3_passive(EtPL,uL),         EtR, SR, cpu_flux3_passive(EtPR,uR));
	pnc4 = cpu_fhll_passive(phicL, SL, cpu_flux4_passive(rhoL,uL,phicL), phicR, SR, cpu_flux4_passive(rhoR,uR,phicR));
	pnc5 = cpu_fhll_passive(phincL, SL, 0, phincR, SR, 0) - 0.5*(SR+SL)/(SR-SL)*vnc5;
      }
      else if (SR < 0){
	pnc1 = cpu_flux1_passive(rhoR,uR);
	pnc2 = cpu_flux2_passive(rhoR,uR,pR);
	pnc3 = cpu_flux3_passive(EtPR,uR);
	pnc4 = cpu_flux4_passive(rhoR,uR,phicR);
	pnc5 = 0.5*vnc5;
      }

      // first
      scalar qL = -pnc1;
      q[(t*N_F+0)*2+0] = qL;
      q[(t*N_F+0)*2+1] = -qL;
      
      //second: fx = rho*u*u+Bx*Bx+Pbar; 
      qL = -pnc2;
      q[(t*N_F+1)*2+0] = qL;
      q[(t*N_F+1)*2+1] = -qL;

      //third: fx = EtplusP*u; 
      qL = -pnc3;
      q[(t*N_F+2)*2+0] = qL; 
      q[(t*N_F+2)*2+1] = -qL;

      //fourth: fx = rho*u*phic
      qL = -pnc4;
      q[(t*N_F+3)*2+0] = qL;
      q[(t*N_F+3)*2+1] = -qL;
      
      //fifth: 
      	qL = -pnc5;
      q[(t*N_F+4)*2+0] = qL - 0.5*vnc5; 
      q[(t*N_F+4)*2+1] = -qL- 0.5*vnc5;
    }

    // Non-conservative Roe flux
    else if (flux == 2){
      printf("Roe flux not implemented\n");
    }
  }
}


void cpu_redistribute_sf(int D, int N_G, int N_E, int N_F, scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac){

  for(int e = 0; e < N_E; e++){
    scalar j = J[e];
    for(int g = 0; g < N_G; g++){
      for(int fc = 0; fc < N_F; fc++){
	scalar sol = 0.0; 

	sJ[(e*N_F+fc)*N_G+g] = s[(e*N_F+fc)*N_G+g] * j;
	for(int alpha = 0; alpha < D; alpha++){
	  for(int a = 0; a < D; a++){
	    sol += invJac[((e*N_G+g)*D+alpha)*D+a]*f[((e*N_F+fc)*N_G+g)*D+a] * j;
	  }
	  fJ[((e*N_F+fc)*N_G+g)*D+alpha] = sol;
	  sol = 0;
	}
      }
    }
  }
}

void cpu_gemm_sf(int D, int N_G, int N_s, int N_E, int N_F, scalar* S, scalar* F, scalar* sJ, scalar* fJ, scalar* phi_w, scalar* dphi_w){

  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      for(int fc = 0; fc < N_F; fc++){
	scalar sol = 0.0;

	// S = phi_w.transpose() x sJ
	for(int g = 0; g < N_G; g++){
	  sol += phi_w[i*N_G+g] * sJ[(e*N_F+fc)*N_G+g];
	}
	S[(e*N_F+fc)*N_s+i] = sol;
	sol = 0.0;
	
	// F = dphi_w.transpose() x fJ
	sol = 0.0; 
	for(int g = 0; g < N_G; g++){
	  for(int a = 0; a < D; a++){
	    sol += dphi_w[(i*N_G+g)*D+a] * fJ[((e*N_F+fc)*N_G+g)*D+a];
	  }
	}
	F[(e*N_F+fc)*N_s+i] = sol;
	sol = 0.0;
      }
    }
  }
}

void cpu_redistribute_q(int M_G, int M_T, int N_F, scalar* qJ, scalar* q){

  for(int t = 0; t < M_T; t++){
    for(int fc = 0; fc < N_F; fc++){
      qJ[(t*N_F+fc)*2+0] = q[(t*N_F+fc)*2+0];
      qJ[(t*N_F+fc)*2+1] = q[(t*N_F+fc)*2+1];
    }
  }
}

 

 void cpu_gemm_q(int M_G, int M_s, int M_T, int N_F, scalar* Qtcj, scalar* qJ, scalar* psi_w){

  for(int t = 0; t < M_T; t++){
    for(int j = 0; j < M_s; j++){
      for(int fc = 0; fc < N_F; fc++){
	scalar sol = 0.0;
	
	// Qtcj = psi_w.transpose() x qJ
	for(int d = 0; d < 2; d++){
	  for(int g = 0; g < M_G; g++){
	    sol += psi_w[j*M_G+g] * qJ[((t*N_F+fc)*2+d)*M_G+g];
	  }
	  Qtcj[((t*N_F+fc)*2+d)*M_s+j] = sol;
	  sol = 0.0;
	}
      }
    }
  }
}

void cpu_solve(int N_s, int N_E, int N_F, scalar* DU, scalar* S, scalar* F, scalar* Q, scalar* Minv, scalar Dt){

  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      for(int fc = 0; fc < N_F; fc++){
	scalar sol = 0.0;
	
	for(int ii = 0; ii < N_s; ii++){
	  sol += Minv[(e*N_s+ii)*N_s+i]*(S[(e*N_F+fc)*N_s+ii] + F[(e*N_F+fc)*N_s+ii] + Q[(e*N_F+fc)*N_s+ii]);
	}
	DU[(e*N_F+fc)*N_s+i] = Dt*sol;
	sol = 0.0;
      }
    }
  }
}

void cpu_average_cell_p0(const int N_s, const int N_E, const int N_F, scalar* DU){

  scalar average = 0.0;
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      average = 0.0;
      for(int i = 0; i < N_s; i++){
	average += DU[(e*N_F+fc)*N_s+i];
      }
      average = average/N_s;
      for(int i = 0; i < N_s; i++){
	DU[(e*N_F+fc)*N_s+i] = average;
      }
    }
  }
}


void cpu_zeroVector(int N_s, int N_E, int N_F, scalar* Q){

  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      for(int fc = 0; fc < N_F; fc++){
	Q[(e*N_F+fc)*N_s+i] = 0.0;
      }
    }
  }
}


//===============================================================
//
//  Host C functions
//
//===============================================================

extern "C" 
void Lcpu_equal(int N_s, int N_E, int N_F, scalar* A, scalar* B){
  cpu_equal(N_s, N_E, N_F, A, B);
}

extern "C" 
void Lcpu_add(int N_s, int N_E, int N_F, scalar* A, scalar* B, scalar c){
  cpu_add(N_s, N_E, N_F, A, B, c);
}

extern "C" 
void Lcpu_mapToFace_shallow(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){
  cpu_mapToFace_shallow(M_s, M_T, N_F, map, U, UF);
}

extern "C" 
void Lcpu_mapToFace_mhd(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){
  cpu_mapToFace_mhd(M_s, M_T, N_F, map, U, UF);
}

extern "C" 
void Lcpu_mapToFace_multifluid(int M_s, int M_T, int N_F, int N_s, int boundaryMap, scalar* U, scalar* UF){
  cpu_mapToFace_multifluid(M_s, M_T, N_F, N_s, boundaryMap, U, UF);
}

extern "C" 
void Lcpu_mapToFace_passive(int M_s, int M_T, int N_F, int N_s, int boundaryMap, scalar* U, scalar* UF){
  cpu_mapToFace_passive(M_s, M_T, N_F, N_s, boundaryMap, U, UF);
}

extern "C"
void Lcpu_boundary(int M_s, int N_F, int M_B, int* boundaryMap, scalar* UF){
  cpu_boundary(M_s, N_F, M_B, boundaryMap, UF);
}

extern "C" 
void Lcpu_mapToElement(int N_s, int N_E, int N_F, scalar* Q, scalar* q){
  cpu_mapToElement(N_s, N_E, N_F, Q, q);
}

extern "C" 
void Lcpu_collocationU(int D, int N_G, int N_s, int N_E, int N_F, scalar* Ug, scalar* dUg, scalar* phi, scalar* dphi, scalar* U){
  cpu_collocationU(D, N_G, N_s, N_E, N_F, Ug, dUg, phi, dphi, U);
}

extern "C" 
void Lcpu_collocationUF(int M_G, int M_s, int M_T, int N_F, scalar* UgF, scalar* psi, scalar* UF){
  cpu_collocationUF(M_G, M_s, M_T, N_F, UgF, psi, UF);
}

extern "C" 
void Lcpu_evaluate_sf_shallow(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0){
  cpu_evaluate_sf_shallow(D, N_G, N_E, N_F, s, f, Ug, H0, G0);
}

extern "C" 
void Lcpu_evaluate_sf_mhd(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma){
  cpu_evaluate_sf_mhd(D, N_G, N_E, N_F, s, f, Ug, dUg, invJac, gamma);
}

extern "C" 
void Lcpu_evaluate_sf_multifluid(int D, int N_G, int N_E, int N_F, int model, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){
  cpu_evaluate_sf_multifluid(D, N_G, N_E, N_F, model, s, f, Ug, dUg, invJac);
}

extern "C" 
void Lcpu_evaluate_sf_passive(int D, int N_G, int N_E, int N_F, scalar gamma, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){
  cpu_evaluate_sf_passive(D, N_G, N_E, N_F, gamma, s, f, Ug, dUg, invJac);
}

extern "C" 
void Lcpu_evaluate_q_shallow(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals){
  cpu_evaluate_q_shallow(M_G, M_T, N_F, q, UgF, H0, G0, normals);
}

extern "C" 
void Lcpu_evaluate_q_mhd(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar gamma, scalar* normals){
  cpu_evaluate_q_mhd(M_G, M_T, N_F, q, UgF, gamma, normals);
}

extern "C" 
void Lcpu_evaluate_q_multifluid(int M_G, int M_T, int N_F, int flux, int model, scalar* q, scalar* UgF){
  cpu_evaluate_q_multifluid(M_G, M_T, N_F, flux, model, q, UgF);
}

extern "C" 
void Lcpu_evaluate_q_passive(int M_G, int M_T, int N_F, int flux, scalar gamma, scalar* q, scalar* UgF){
  cpu_evaluate_q_passive(M_G, M_T, N_F, flux, gamma, q, UgF);
}

extern "C" 
void Lcpu_redistribute_sf(int D, int N_G, int N_E, int N_F, scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac){
  cpu_redistribute_sf(D, N_G, N_E, N_F, sJ, fJ, s, f, J, invJac);
}

extern "C" 
void Lcpu_gemm_sf(int D, int N_G, int N_s, int N_E, int N_F, scalar* S, scalar* F, scalar* sJ, scalar* fJ, scalar* phi_w, scalar* dphi_w){
  cpu_gemm_sf(D, N_G, N_s, N_E, N_F, S, F, sJ, fJ, phi_w, dphi_w);
}

extern "C"
void Lcpu_redistribute_q(int M_G, int M_T, int N_F, scalar* qJ, scalar* q){
  cpu_redistribute_q(M_G, M_T, N_F, qJ, q);
}

extern "C" 
void Lcpu_gemm_q(int M_G, int M_s, int M_T, int N_F, scalar* Qtcj, scalar* qJ, scalar* psi_w){
  cpu_gemm_q(M_G, M_s, M_T, N_F, Qtcj, qJ, psi_w);
}

extern "C" 
void Lcpu_solve(int N_s, int N_E, int N_F, scalar* DU, scalar* S, scalar* F, scalar* Q, scalar* Minv, scalar Dt){
  cpu_solve(N_s, N_E, N_F, DU, S, F, Q, Minv, Dt);
}

extern "C"
void Lcpu_average_cell_p0(const int N_s, const int N_E, const int N_F, scalar* DU){
  cpu_average_cell_p0(N_s, N_E, N_F, DU);
}
