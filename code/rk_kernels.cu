#include <rk_kernels.h>
#include <cstdlib>
#include <stdio.h>
#include <upa.h>
#include <constants.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void solve(int N_s, int N_E, scalar Dt, scalar* Minv, scalar* f, scalar* DU){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
    int fc= threadIdx.y;
#endif

  scalar sol = 0.0;
	
  for(int ii = 0; ii < N_s; ii++){
    sol += Minv[(e*N_s+ii)*N_s+i]*f[(e*N_F+fc)*N_s+ii];
  }
  DU[(e*N_F+fc)*N_s+i] = Dt*sol;
  sol = 0.0;

#ifdef USE_CPU
      }
    }
#endif
  }
}


//==========================================================================
arch_global void average_cell_p0(const int N_s, const int N_E,  scalar* DU){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int fc= threadIdx.y;
#endif
  
  scalar average = 0.0;
  for(int i = 0; i < N_s; i++){
    average += DU[(e*N_F+fc)*N_s+i];
  }
  average = average/N_s;
  for(int i = 0; i < N_s; i++){
    DU[(e*N_F+fc)*N_s+i] = average;
  }

#ifdef USE_CPU
    }
#endif
  }
}

//==========================================================================
arch_global void findUPA(const int N_s, const int N_E,  scalar* U, scalar* UPA){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
#endif

#ifdef PASSIVE
#ifdef ONED
  UPA[e*N_s+i] = oned_passive_upa(U[(e*N_F+0)*N_s+i],                    // rho
				  U[(e*N_F+1)*N_s+i]/U[(e*N_F+0)*N_s+i], // u
				  U[(e*N_F+2)*N_s+i],                    // E
				  constants::GLOBAL_GAMMA);              // gamma
#elif TWOD
  UPA[e*N_s+i] = twod_passive_upa(U[(e*N_F+0)*N_s+i],                    // rho
				  U[(e*N_F+1)*N_s+i]/U[(e*N_F+0)*N_s+i], // u
				  U[(e*N_F+2)*N_s+i]/U[(e*N_F+0)*N_s+i], // v
				  U[(e*N_F+3)*N_s+i],                    // E
				  constants::GLOBAL_GAMMA);              // gamma
#endif // dimensions
#elif MULTIFLUID
#ifdef ONED
  UPA[e*N_s+i] = oned_multifluid_upa(U[(e*N_F+0)*N_s+i],                    // rho
				     U[(e*N_F+1)*N_s+i]/U[(e*N_F+0)*N_s+i], // u
				     U[(e*N_F+2)*N_s+i],                    // E
				     U[(e*N_F+3)*N_s+i]);                   // alpha
#elif TWOD
  UPA[e*N_s+i] = twod_multifluid_upa(U[(e*N_F+0)*N_s+i],                    // rho
				     U[(e*N_F+1)*N_s+i]/U[(e*N_F+0)*N_s+i], // u
				     U[(e*N_F+2)*N_s+i]/U[(e*N_F+0)*N_s+i], // v
				     U[(e*N_F+3)*N_s+i],                    // E
				     U[(e*N_F+4)*N_s+i]);                   // alpha
#endif //dimesions
#elif STIFFENED
#ifdef ONED
  UPA[e*N_s+i] = oned_stiffened_upa(U[(e*N_F+0)*N_s+i],                    // rho
				    U[(e*N_F+1)*N_s+i]/U[(e*N_F+0)*N_s+i], // u
				    U[(e*N_F+2)*N_s+i],                    // E
				    U[(e*N_F+3)*N_s+i],                    // alpha
  				    U[(e*N_F+4)*N_s+i]);                   // beta
#elif TWOD
  UPA[e*N_s+i] = twod_stiffened_upa(U[(e*N_F+0)*N_s+i],                    // rho
				    U[(e*N_F+1)*N_s+i]/U[(e*N_F+0)*N_s+i], // u
				    U[(e*N_F+2)*N_s+i]/U[(e*N_F+0)*N_s+i], // v
				    U[(e*N_F+3)*N_s+i],                    // E
				    U[(e*N_F+4)*N_s+i],                    // alpha
  				    U[(e*N_F+5)*N_s+i]);                   // beta
  
#endif //dimesions
#endif // problem type

#ifdef USE_CPU
    }
#endif
  }
}

//==========================================================================
//
//  Host C functions
//
//==========================================================================

extern "C" 
void Lsolver(int N_s, int N_E, scalar Dt, scalar* Minv, scalar* f, scalar* DU){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  solve arch_args (N_s, N_E, Dt, Minv, f, DU);
};

extern "C"
void Laverage_cell_p0(const int N_s, const int N_E,  scalar* DU){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(1,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  average_cell_p0 arch_args (N_s, N_E, DU);
}

extern "C"
void LfindUPA(const int N_s, const int N_E, scalar* U, scalar* UPA){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif
    
    findUPA arch_args (N_s, N_E, U, UPA);
}
