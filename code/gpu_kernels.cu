#include <gpu_kernels.h>
// Kernel definitions
__global__ void gpu_equal(int N_s, int N_E, int N_F, scalar* A, scalar* B){

  int e = blockIdx.x;
  int i = threadIdx.x;
  int fc = threadIdx.y;

  A[(e*N_F+fc)*N_s+i] = B[(e*N_F+fc)*N_s+i];
}

__global__ void gpu_add(int N_s, int N_E, int N_F, scalar* A, scalar* B, scalar c){

  // A = A + c*B
  int e = blockIdx.x;
  int i = threadIdx.x;
  int fc = threadIdx.y;

  A[(e*N_F+fc)*N_s+i] =  A[(e*N_F+fc)*N_s+i] + c*B[(e*N_F+fc)*N_s+i];
}


__global__ void gpu_mapToFace_shallow(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){

  int t = blockIdx.x;
  int j = threadIdx.x;
  int fc= threadIdx.y;
  
  int idx = -1;
  int face;
  
  for(int d = 0; d < 2; d++){
    face= ((t*N_F+fc)*2+d)*M_s+j;
    idx = map[face];
    if(idx != -1){
      UF[face] = U[idx];
    }
    else if (idx == -1){
      if      (fc == 0) UF[((t*N_F+fc)*2+1)*M_s+j] = UF[((t*N_F+fc)*2+0)*M_s+j]; // eta
      else if (fc == 1) UF[((t*N_F+fc)*2+1)*M_s+j] =-UF[((t*N_F+fc)*2+0)*M_s+j]; // ux
      else if (fc == 2) UF[((t*N_F+fc)*2+1)*M_s+j] =-UF[((t*N_F+fc)*2+0)*M_s+j]; // uy
    }
  }
}

__global__ void gpu_mapToFace_mhd(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){

  int t = blockIdx.x;
  int j = threadIdx.x;
  int fc= threadIdx.y;
  
  int idx = -1;
  int face;
  
  for(int d = 0; d < 2; d++){
    face= ((t*N_F+fc)*2+d)*M_s+j;
    idx = map[face];
    if(idx != -1){
      UF[face] = U[idx];
    }
  }
}

__global__ void gpu_mapToFace_multifluid(int M_s, int M_T, int N_F, int N_s, int boundaryMap, scalar* U, scalar* UF){

  int t = blockIdx.x;
  int fc= threadIdx.y;
    
  // Start boundary
  if (t==0){
    UF[(0*N_F+fc)*2+1]       = U[(0*N_F+fc)*N_s+0];         
    if      (boundaryMap == 0){      //farfield
      UF[(0*N_F+fc)*2+0]        = UF[(0*N_F+fc)*2+1];}
    else if (boundaryMap == M_T-1){  //periodic
      UF[(0*N_F+fc)*2+0]        = UF[((M_T-1)*N_F+fc)*2+0];}
  }
  // End boundary
  else if (t==M_T-1){
    UF[((M_T-1)*N_F+fc)*2+0] = U[((M_T-2)*N_F+fc)*N_s+1];
    if      (boundaryMap == 0){      //farfield
      UF[((M_T-1)*N_F+fc)*2+1]  = UF[((M_T-1)*N_F+fc)*2+0];}
    else if (boundaryMap == M_T-1){  //periodic
      UF[((M_T-1)*N_F+fc)*2+1]  = UF[(0*N_F+fc)*2+1];}
  }
  // All other boundaries
  else{
    UF[(t*N_F+fc)*2+0] = U[((t-1)*N_F+fc)*N_s+1];
    UF[(t*N_F+fc)*2+1] = U[(t*N_F+fc)*N_s+0];
  }
}

__global__ void gpu_boundary(int M_s, int N_F, int M_B, int* boundaryMap, scalar* UF){

  int t1 = boundaryMap[blockIdx.x*2+0];
  int t2 = boundaryMap[blockIdx.x*2+1];
  int j = threadIdx.x;
  int fc= threadIdx.y;
  UF[((t1*N_F+fc)*2+1)*M_s+j] = UF[((t2*N_F+fc)*2+0)*M_s+j]; 
}

__global__ void gpu_mapToElement(int N_s, int N_E, int N_F, scalar* Q, scalar* q){

  int e = blockIdx.x;
  int fc= threadIdx.y;

  Q[(e*N_F+fc)*N_s+0] = q[(e*N_F+fc)*2+1];
  Q[(e*N_F+fc)*N_s+1] = q[((e+1)*N_F+fc)*2+0];
}

__global__ void gpu_collocationU(int D, int N_G, int N_s, int N_E, int N_F, scalar* Ug, scalar* dUg, scalar* phi, scalar* dphi, scalar* U){

  int e = blockIdx.x;
  int g = threadIdx.x;
  int fc= threadIdx.y;

  scalar sol = 0;

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

__global__ void gpu_collocationUF(int M_G, int M_s, int M_T, int N_F, scalar* UgF, scalar* psi, scalar* UF){

  int t = blockIdx.x;
  int g = threadIdx.x;
  int fc= threadIdx.y;

  scalar sol = 0;
  for(int d = 0; d < 2; d++){
    for(int j = 0; j < M_s; j++){
      sol += psi[j*M_G+g] * UF[((t*N_F+fc)*2+d)*M_s+j];
    }
    UgF[((t*N_F+fc)*2+d)*M_G+g] = sol;
    sol = 0.0;
  }
}

__global__ void gpu_evaluate_sf_shallow(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0){

  int e = blockIdx.x;
  int g = threadIdx.x;

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

__device__ scalar gpu_flux1_mhd(scalar rho, scalar u){return rho*u;}                                        // for f0X, f0X, f0X
__device__ scalar gpu_flux2_mhd(scalar rho, scalar u, scalar Bx, scalar pbar){return rho*u*u-Bx*Bx+pbar;}   // for f1X, f2Y, f3Z
__device__ scalar gpu_flux3_mhd(scalar rho, scalar u, scalar v, scalar Bx, scalar By){return rho*u*v-Bx*By;}// for f1Y, f1Z, f2X, f2Z, f3X, f3Y
__device__ scalar gpu_flux4_mhd(scalar u, scalar By, scalar v, scalar Bx){return u*By-v*Bx;}                // for f4Y, f4Z, f5X, f5Z, f6X, f6Y
__device__ scalar gpu_flux5_mhd(scalar EtplusPbar, scalar u, scalar vdotB, scalar Bx) {return EtplusPbar*u - vdotB*Bx;} // for f7X, f7Y, f7Z


__global__ void gpu_evaluate_sf_mhd(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma){

  int e = blockIdx.x;
  int g = threadIdx.x;

  scalar rho = Ug[(e*N_F+0)*N_G+g];   
  scalar u   = Ug[(e*N_F+1)*N_G+g]/rho; 
  scalar v   = Ug[(e*N_F+2)*N_G+g]/rho; 
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
  f[((e*N_F+0)*N_G+g)*D+0] = gpu_flux1_mhd(rho,u);                //rho*u; 
  f[((e*N_F+1)*N_G+g)*D+0] = gpu_flux2_mhd(rho,u,Bx,Pbar);        //rho*u*u-Bx*Bx+Pbar; 
  f[((e*N_F+2)*N_G+g)*D+0] = gpu_flux3_mhd(rho,u,v,Bx,By);        //rho*u*v-Bx*By; 
  f[((e*N_F+3)*N_G+g)*D+0] = 0;                                   //0;
  f[((e*N_F+4)*N_G+g)*D+0] = gpu_flux4_mhd(u,By,v,Bx);            //u*By-v*Bx;
  f[((e*N_F+5)*N_G+g)*D+0] = gpu_flux5_mhd(EtPbar,u,vdotB,Bx);    //EtplusPbar*u-vdotB*Bx;
      
  // Flux derive par rapport a y
  f[((e*N_F+0)*N_G+g)*D+1] = gpu_flux1_mhd(rho,v);                //rho*v;
  f[((e*N_F+1)*N_G+g)*D+1] = gpu_flux3_mhd(rho,v,u,By,Bx);        //rho*v*u-By*Bx;
  f[((e*N_F+2)*N_G+g)*D+1] = gpu_flux2_mhd(rho,v,By,Pbar);        //rho*v*v-By*By+Pbar;
  f[((e*N_F+3)*N_G+g)*D+1] = gpu_flux4_mhd(v,Bx,u,By);            //v*Bx-u*By;
  f[((e*N_F+4)*N_G+g)*D+1] = 0;                                   //0;
  f[((e*N_F+5)*N_G+g)*D+1] = gpu_flux5_mhd(EtPbar,v,vdotB,By);    //EtplusPbar*v-vdotB*By;

}

__device__ scalar gpu_flux1_multifluid(scalar rho, scalar u){return rho*u;}                  
__device__ scalar gpu_flux2_multifluid(scalar rho, scalar u, scalar p){return rho*u*u+p;} 
__device__ scalar gpu_flux3_multifluid(scalar EtplusP, scalar u) {return EtplusP*u;}
__device__ scalar gpu_flux4_multifluid(scalar rho, scalar u, scalar gamma) {return rho*u*gamma;}
 
__global__ void gpu_evaluate_sf_multifluid(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){

  int e = blockIdx.x;
  int g = threadIdx.x;

  s[(e*N_F+0)*N_G+g] = 0;
  s[(e*N_F+1)*N_G+g] = 0;
  s[(e*N_F+2)*N_G+g] = 0;
  s[(e*N_F+3)*N_G+g] = 0;

  scalar rho = Ug[(e*N_F+0)*N_G+g];   
  scalar u   = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
  scalar Et  = Ug[(e*N_F+2)*N_G+g];
  scalar gamma = Ug[(e*N_F+3)*N_G+g]/rho;
  scalar p = (gamma-1)*(Et - 0.5*rho*u*u);
  scalar EtplusP = Et + p;

  // Flux derive par rapport a x
  f[((e*N_F+0)*N_G+g)*D+0] = gpu_flux1_multifluid(rho,u);        //rho*u; 
  f[((e*N_F+1)*N_G+g)*D+0] = gpu_flux2_multifluid(rho,u,p);      //rho*u*u+P; 
  f[((e*N_F+2)*N_G+g)*D+0] = gpu_flux3_multifluid(EtplusP,u);    //EtplusP*u;
  f[((e*N_F+3)*N_G+g)*D+0] = gpu_flux4_multifluid(rho,u,gamma); 
}


__global__ void gpu_evaluate_q_shallow(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals){

  int t = blockIdx.x;
  int g = threadIdx.x;

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

__global__ void gpu_evaluate_q_mhd(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar gamma, scalar* normals){

  int t = blockIdx.x;
  int g = threadIdx.x;

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
  extern __shared__ scalar vap[];
  scalar alfenL = BdotnL/sqrt(rhoL);
  scalar alfenR = BdotnR/sqrt(rhoR);
  scalar a2L = (gamma*pL+BdotBL)/rhoL;
  scalar a2R = (gamma*pR+BdotBR)/rhoR;
  scalar cfL = sqrt(0.5*(a2L + sqrt( a2L*a2L - 4.0*gamma*pL*BdotnL*BdotnL/(rhoL*rhoL))));
  scalar cfR = sqrt(0.5*(a2R + sqrt( a2R*a2R - 4.0*gamma*pR*BdotnR*BdotnR/(rhoR*rhoR))));
  scalar csL = sqrt(0.5*(a2L - sqrt( a2L*a2L - 4.0*gamma*pL*BdotnL*BdotnL/(rhoL*rhoL))));
  scalar csR = sqrt(0.5*(a2R - sqrt( a2R*a2R - 4.0*gamma*pR*BdotnR*BdotnR/(rhoR*rhoR))));
      
  vap[(g*2+0)*sizevap+0] = fabs(vdotnL);
  vap[(g*2+0)*sizevap+1] = fabs(vdotnL) + alfenL;
  vap[(g*2+0)*sizevap+2] = fabs(vdotnL) - alfenL;
  vap[(g*2+0)*sizevap+3] = fabs(vdotnL) + cfL;
  vap[(g*2+0)*sizevap+4] = fabs(vdotnL) - cfL;
  vap[(g*2+0)*sizevap+5] = fabs(vdotnL) + csL;
  vap[(g*2+0)*sizevap+6] = fabs(vdotnL) - csL;
  vap[(g*2+0)*sizevap+7] = fabs(vdotnL);
      
  vap[(g*2+1)*sizevap+0] = fabs(vdotnR);
  vap[(g*2+1)*sizevap+1] = fabs(vdotnR) + alfenR;
  vap[(g*2+1)*sizevap+2] = fabs(vdotnR) - alfenR;
  vap[(g*2+1)*sizevap+3] = fabs(vdotnR) + cfR;
  vap[(g*2+1)*sizevap+4] = fabs(vdotnR) - cfR;
  vap[(g*2+1)*sizevap+5] = fabs(vdotnR) + csR;
  vap[(g*2+1)*sizevap+6] = fabs(vdotnR) - csR;
  vap[(g*2+1)*sizevap+7] = fabs(vdotnR);

  scalar maxvap = 0;
  for (int k = 0; k < 2*sizevap; k++){
    if (maxvap<vap[g*16+k]) maxvap = vap[g*16+k];
  }


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
  scalar qL = -0.5*((gpu_flux1_mhd(rhoL,uL) + gpu_flux1_mhd(rhoR,uR))*nx +
		    (gpu_flux1_mhd(rhoL,vL) + gpu_flux1_mhd(rhoR,vR))*ny 
		    -maxvap*(rhoR-rhoL));
  q[((t*N_F+0)*2+0)*M_G+g] = qL;
  q[((t*N_F+0)*2+1)*M_G+g] = -qL;
      
  //second: fx = rho*u*u+Bx*Bx+Pbar; fy = rho*v*u-By*Bx; fz = rho*w*u-Bz*Bx;
  qL = -0.5*((gpu_flux2_mhd(rhoL,uL,BxL,PbarL)  + gpu_flux2_mhd(rhoR,uR,BxR,PbarR) )*nx +
	     (gpu_flux3_mhd(rhoL,vL,uL,ByL,BxL) + gpu_flux3_mhd(rhoR,vR,uR,ByR,BxR))*ny  
	     -maxvap*(rhoR*uR-rhoL*uL));
  q[((t*N_F+1)*2+0)*M_G+g] = qL  - BxL*upBdotnL;
  q[((t*N_F+1)*2+1)*M_G+g] = -qL + BxR*upBdotnR;

  //third: fx = rho*u*v-Bx*By; fy = rho*v*v-By*By+Pbar; fz = rho*w*v-Bz*By;
  qL = -0.5*((gpu_flux3_mhd(rhoL,uL,vL,BxL,ByL) + gpu_flux3_mhd(rhoR,uR,vR,BxR,ByR))*nx +
	     (gpu_flux2_mhd(rhoL,vL,ByL,PbarL)  + gpu_flux2_mhd(rhoR,vR,ByR,PbarR) )*ny  
	     -maxvap*(rhoR*vR-rhoL*vL));
  q[((t*N_F+2)*2+0)*M_G+g] = qL  - ByL*upBdotnL;
  q[((t*N_F+2)*2+1)*M_G+g] = -qL + ByR*upBdotnR;
      
  //fourth: fx = 0;  fy = v*Bx-u*By; fz = w*Bx-u*Bz;
  qL = -0.5*((0                            + 0                           )*nx +
	     (gpu_flux4_mhd(vL,BxL,uL,ByL) + gpu_flux4_mhd(vR,BxR,uR,ByR))*ny 
	     -maxvap*(BxR-BxL));
  q[((t*N_F+3)*2+0)*M_G+g] = qL  - uL*upBdotnL;
  q[((t*N_F+3)*2+1)*M_G+g] = -qL + uR*upBdotnR;

  //fifth: fx = u*By-v*Bx;  fy = 0; fz = w*By-v*Bz;
  qL = -0.5*((gpu_flux4_mhd(uL,ByL,vL,BxL) + gpu_flux4_mhd(uR,ByR,vR,BxR))*nx +
	     (0                            + 0                           )*ny 
	     -maxvap*(ByR-ByL));
  q[((t*N_F+4)*2+0)*M_G+g] = qL  - vL*upBdotnL;
  q[((t*N_F+4)*2+1)*M_G+g] = -qL + vR*upBdotnR;

  //sixth: fx = EtplusPbar*u-vdotB*Bx; fy = EtplusPbar*v-vdotB*By; fz = EtplusPbar*w-vdotB*Bz;
  qL = -0.5*((gpu_flux5_mhd(EtPbarL,uL,vdotBL,BxL) + gpu_flux5_mhd(EtPbarR,uR,vdotBR,BxR))*nx +
	     (gpu_flux5_mhd(EtPbarL,vL,vdotBL,ByL) + gpu_flux5_mhd(EtPbarR,vR,vdotBR,ByR))*ny 
	     -maxvap*(EtR-EtL));
  q[((t*N_F+5)*2+0)*M_G+g] = qL  - vdotBL*upBdotnL;
  q[((t*N_F+5)*2+1)*M_G+g] = -qL + vdotBR*upBdotnR; 

}


__global__ void gpu_evaluate_q_multifluid(int M_G, int M_T, int N_F, scalar* q, scalar* UgF){

  int t = blockIdx.x;

  scalar rhoL= UgF[(t*N_F+0)*2+0];
  scalar rhoR= UgF[(t*N_F+0)*2+1];
  scalar uL  = UgF[(t*N_F+1)*2+0]/rhoL;
  scalar uR  = UgF[(t*N_F+1)*2+1]/rhoR;
  scalar EtL = UgF[(t*N_F+2)*2+0];
  scalar EtR = UgF[(t*N_F+2)*2+1];
  scalar gammaL = UgF[(t*N_F+3)*2+0]/rhoL;
  scalar gammaR = UgF[(t*N_F+3)*2+1]/rhoR;
  scalar pL = (gammaL-1)*(EtL - 0.5*rhoL*uL*uL);
  scalar pR = (gammaR-1)*(EtR - 0.5*rhoR*uR*uR);
  scalar EtPL = EtL+pL;
  scalar EtPR = EtR+pR;

  // Evaluate the right and left eigenvalues
  int sizevap = 4;
  extern __shared__ scalar vap[];
  scalar aL = sqrt((gammaL*pL)/rhoL);
  scalar aR = sqrt((gammaR*pR)/rhoR);

  vap[0*sizevap+0] = fabs(uL);
  vap[0*sizevap+1] = fabs(uL) + aL;
  vap[0*sizevap+2] = fabs(uL) - aL;
  vap[0*sizevap+3] = fabs(uL);
      
  vap[1*sizevap+0] = fabs(uR);
  vap[1*sizevap+1] = fabs(uR) + aR;
  vap[1*sizevap+2] = fabs(uR) - aR;
  vap[1*sizevap+3] = fabs(uR);

  scalar maxvap = 0;
  for (int k = 0; k < 2*sizevap; k++){
    if (maxvap<vap[k]) maxvap = vap[k];
  }
      
  //
  // Evaluate the fluxes on the right and left
  //

  //first: fx = rho*u; 
  scalar qL = -0.5*(gpu_flux1_multifluid(rhoL,uL) + gpu_flux1_multifluid(rhoR,uR)
		    -maxvap*(rhoR-rhoL));
  q[(t*N_F+0)*2+0] = qL;
  q[(t*N_F+0)*2+1] = -qL;
      
  //second: fx = rho*u*u+Bx*Bx+Pbar; 
  qL = -0.5*(gpu_flux2_multifluid(rhoL,uL,pL)  + gpu_flux2_multifluid(rhoR,uR,pR)
	     -maxvap*(rhoR*uR-rhoL*uL));
  q[(t*N_F+1)*2+0] = qL;
  q[(t*N_F+1)*2+1] = -qL;

  //third: fx = EtplusP*u; 
  qL = -0.5*(gpu_flux3_multifluid(EtPL,uL) + gpu_flux3_multifluid(EtPR,uR)
	     -maxvap*(EtR-EtL));
  q[(t*N_F+2)*2+0] = qL; 
  q[(t*N_F+2)*2+1] = -qL;

  //fourth: fx = rho*u*gamma; 
  qL = -0.5*(gpu_flux4_multifluid(rhoL,uL,gammaL) + gpu_flux4_multifluid(rhoR,uR,gammaR)
	     -maxvap*(rhoR*gammaR-rhoL*gammaL));
  q[(t*N_F+3)*2+0] = qL; 
  q[(t*N_F+3)*2+1] = -qL;
  
}


__global__ void gpu_redistribute_sf(int D, int N_G, int N_E, int N_F, scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac){

  int e = blockIdx.x;
  int g = threadIdx.x;
  int fc= threadIdx.y;

  scalar j = J[e];
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

__global__ void gpu_gemm_sf(int D, int N_G, int N_s, int N_E, int N_F, scalar* S, scalar* F, scalar* sJ, scalar* fJ, scalar* phi_w, scalar* dphi_w){

  int e = blockIdx.x;
  int i = threadIdx.x;
  int fc= threadIdx.y;

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

__global__ void gpu_redistribute_q(int M_G, int M_T, int N_F, scalar* qJ, scalar* q, scalar* JF){

  int t = blockIdx.x;
  int g = threadIdx.x;
  int fc= threadIdx.y;

  for(int d = 0; d < 2; d++){
    qJ[((t*N_F+fc)*2+d)*M_G+g] = q[((t*N_F+fc)*2+d)*M_G+g] * JF[t*2+d];
  }
}

 

__global__  void gpu_gemm_q(int M_G, int M_s, int M_T, int N_F, scalar* Qtcj, scalar* qJ, scalar* psi_w){

  int t = blockIdx.x;
  int j = threadIdx.x;
  int fc= threadIdx.y;

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

__global__ void gpu_solve(int N_s, int N_E, int N_F, scalar* DU, scalar* S, scalar* F, scalar* Q, scalar* Minv, scalar Dt){

  int e = blockIdx.x;
  int i = threadIdx.x;
  int fc= threadIdx.y;

  scalar sol = 0.0;

  for(int ii = 0; ii < N_s; ii++){
    sol += Minv[(e*N_s+ii)*N_s+i]*(S[(e*N_F+fc)*N_s+ii] + F[(e*N_F+fc)*N_s+ii] + Q[(e*N_F+fc)*N_s+ii]);
  }
  DU[(e*N_F+fc)*N_s+i] = Dt*sol;
  sol = 0.0;
}

__global__ void gpu_average_cell_p0(const int N_s, const int N_E, const int N_F, scalar* DU){

  int e = blockIdx.x;
  int fc= threadIdx.y;
  scalar average = 0.0;

  for(int i = 0; i < N_s; i++){
    average += DU[(e*N_F+fc)*N_s+i];
  }
  average = average/N_s;
  for(int i = 0; i < N_s; i++){
    DU[(e*N_F+fc)*N_s+i] = average;
  }
}


__global__ void gpu_zeroVector(int N_s, int N_E, int N_F, scalar* Q){

  int e = blockIdx.x;
  int i = threadIdx.x;
  int fc= threadIdx.y;

  Q[(e*N_F+fc)*N_s+i] = 0.0;
}




//===============================================================
//
//  Host C functions
//
//===============================================================

extern "C" 
void Lgpu_equal(int N_s, int N_E, int N_F, scalar* A, scalar* B){

  dim3 dimBlock(N_s,N_F,1);
  dim3 dimGrid(N_E,1);
  gpu_equal<<<dimGrid,dimBlock>>>(N_s, N_E, N_F, A, B);

}

extern "C" 
void Lgpu_add(int N_s, int N_E, int N_F, scalar* A, scalar* B, scalar c){
  
  dim3 dimBlock(N_s,N_F,1);
  dim3 dimGrid(N_E,1);
  gpu_add<<<dimGrid,dimBlock>>>(N_s, N_E, N_F, A, B, c);
}

extern "C" 
void Lgpu_mapToFace_shallow(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){

  dim3 dimBlock(M_s,N_F,1);
  dim3 dimGrid(M_T,1);
  gpu_mapToFace_shallow<<<dimGrid,dimBlock>>>(M_s, M_T, N_F, map, U, UF);
}

extern "C" 
void Lgpu_mapToFace_mhd(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){

  dim3 dimBlock(M_s,N_F,1);
  dim3 dimGrid(M_T,1);
  gpu_mapToFace_mhd<<<dimGrid,dimBlock>>>(M_s, M_T, N_F, map, U, UF);
}

extern "C" 
void Lgpu_mapToFace_multifluid(int M_s, int M_T, int N_F, int N_s, int boundaryMap, scalar* U, scalar* UF){

  dim3 dimBlock(1,N_F,1);
  dim3 dimGrid(M_T,1);
  gpu_mapToFace_multifluid<<<dimGrid,dimBlock>>>(M_s, M_T, N_F, N_s, boundaryMap, U, UF);
}

extern "C" 
void Lgpu_boundary(int M_s, int N_F, int M_B, int* boundaryMap, scalar* UF){

  dim3 dimBlock(M_s,N_F,1);
  dim3 dimGrid(M_B,1);
  gpu_boundary<<<dimGrid,dimBlock>>>(M_s, N_F, M_B, boundaryMap, UF);
}


extern "C" 
void Lgpu_mapToElement(int N_s, int N_E, int N_F, scalar* Q, scalar* q){

  dim3 dimBlock(1,N_F,1);
  dim3 dimGrid(N_E,1);
  gpu_mapToElement<<<dimGrid,dimBlock>>>(N_s, N_E, N_F, Q, q);
}

extern "C" 
void Lgpu_collocationU(int D, int N_G, int N_s, int N_E, int N_F, scalar* Ug, scalar* dUg, scalar* phi, scalar* dphi, scalar* U){

  dim3 dimBlock(N_G,N_F,1);
  dim3 dimGrid(N_E,1);
  gpu_collocationU<<<dimGrid,dimBlock>>>(D, N_G, N_s, N_E, N_F, Ug, dUg, phi, dphi, U);
}

extern "C" 
void Lgpu_collocationUF(int M_G, int M_s, int M_T, int N_F, scalar* UgF, scalar* psi, scalar* UF){
  
  dim3 dimBlock(M_G,N_F,1);
  dim3 dimGrid(M_T,1);
  gpu_collocationUF<<<dimGrid,dimBlock>>>(M_G, M_s, M_T, N_F, UgF, psi, UF);
}

extern "C" 
void Lgpu_evaluate_sf_shallow(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0){

  dim3 dimBlock(N_G,1,1);
  dim3 dimGrid(N_E,1);
  gpu_evaluate_sf_shallow<<<dimGrid,dimBlock>>>(D, N_G, N_E, N_F, s, f, Ug, H0, G0);
}

extern "C" 
void Lgpu_evaluate_sf_mhd(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma){

  dim3 dimBlock(N_G,1,1);
  dim3 dimGrid(N_E,1);
  gpu_evaluate_sf_mhd<<<dimGrid,dimBlock>>>(D, N_G, N_E, N_F, s, f, Ug, dUg, invJac, gamma);
}

extern "C" 
void Lgpu_evaluate_sf_multifluid(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){

  dim3 dimBlock(N_G,1,1);
  dim3 dimGrid(N_E,1);
  gpu_evaluate_sf_multifluid<<<dimGrid,dimBlock>>>(D, N_G, N_E, N_F, s, f, Ug, dUg, invJac);
}

extern "C" 
void Lgpu_evaluate_q_shallow(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals){

  dim3 dimBlock(M_G,1,1);
  dim3 dimGrid(M_T,1);
  gpu_evaluate_q_shallow<<<dimGrid,dimBlock>>>(M_G, M_T, N_F, q, UgF, H0, G0, normals);
}

extern "C" 
void Lgpu_evaluate_q_mhd(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar gamma, scalar* normals){

  dim3 dimBlock(M_G,1,1);
  dim3 dimGrid(M_T,1);
  gpu_evaluate_q_mhd<<<dimGrid,dimBlock,M_G*2*8*sizeof(scalar)>>>(M_G, M_T, N_F, q, UgF, gamma, normals);
}

extern "C" 
void Lgpu_evaluate_q_multifluid(int M_G, int M_T, int N_F, scalar* q, scalar* UgF){

  dim3 dimBlock(1,1,1);
  dim3 dimGrid(M_T,1);
  gpu_evaluate_q_multifluid<<<dimGrid,dimBlock,2*4*sizeof(scalar)>>>(M_G, M_T, N_F, q, UgF);
}

extern "C" 
void Lgpu_redistribute_sf(int D, int N_G, int N_E, int N_F, scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac){

  dim3 dimBlock(N_G,N_F,1);
  dim3 dimGrid(N_E,1);
  gpu_redistribute_sf<<<dimGrid,dimBlock>>>(D, N_G, N_E, N_F, sJ, fJ, s, f, J, invJac);
}

extern "C" 
void Lgpu_gemm_sf(int D, int N_G, int N_s, int N_E, int N_F, scalar* S, scalar* F, scalar* sJ, scalar* fJ, scalar* phi_w, scalar* dphi_w){

  dim3 dimBlock(N_s,N_F,1);
  dim3 dimGrid(N_E,1);
  gpu_gemm_sf<<<dimGrid,dimBlock>>>(D, N_G, N_s, N_E, N_F, S, F, sJ, fJ, phi_w, dphi_w);
}

extern "C"
void Lgpu_redistribute_q(int M_G, int M_T, int N_F, scalar* qJ, scalar* q, scalar* JF){

  dim3 dimBlock(M_G,N_F,1);
  dim3 dimGrid(M_T,1);
  gpu_redistribute_q<<<dimGrid,dimBlock>>>(M_G, M_T, N_F, qJ, q, JF);
}

extern "C" 
void Lgpu_gemm_q(int M_G, int M_s, int M_T, int N_F, scalar* Qtcj, scalar* qJ, scalar* psi_w){

  dim3 dimBlock(M_s,N_F,1);
  dim3 dimGrid(M_T,1);
  gpu_gemm_q<<<dimGrid,dimBlock>>>(M_G, M_s, M_T, N_F, Qtcj, qJ, psi_w);
}

extern "C" 
void Lgpu_solve(int N_s, int N_E, int N_F, scalar* DU, scalar* S, scalar* F, scalar* Q, scalar* Minv, scalar Dt){

 dim3 dimBlock(N_s,N_F,1);
 dim3 dimGrid(N_E,1);
 gpu_solve<<<dimGrid,dimBlock>>>(N_s, N_E, N_F, DU, S, F, Q, Minv, Dt);
}

extern "C"
void Lgpu_average_cell_p0(const int N_s, const int N_E, const int N_F, scalar* DU){

 dim3 dimBlock(1,N_F,1);
 dim3 dimGrid(N_E,1);
 gpu_average_cell_p0<<<dimGrid,dimBlock>>>(N_s, N_E, N_F, DU);
}

