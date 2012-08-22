#include <rk_kernels.h>
#include <cstdlib>
#include <stdio.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void solve(int N_s, int N_E, int N_F, scalar Dt, scalar* Minv, scalar* f, scalar* DU){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int e = blockIdx.x;
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
  }
#endif
}


//==========================================================================
arch_global void average_cell_p0(const int N_s, const int N_E, const int N_F, scalar* DU){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int e = blockIdx.x;
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
  }
#endif
}

//==========================================================================
//
//  Host C functions
//
//==========================================================================

extern "C" 
void Lsolver(int N_s, int N_E, int N_F, scalar Dt, scalar* Minv, scalar* f, scalar* DU){
#ifdef USE_GPU
  dim3 dimBlock(N_s,N_F,1);
  dim3 dimGrid(N_E,1);
#endif

  solve arch_args (N_s, N_E, N_F, Dt, Minv, f, DU);
};

extern "C"
void Laverage_cell_p0(const int N_s, const int N_E, const int N_F, scalar* DU){

#ifdef USE_GPU
  dim3 dimBlock(1,N_F,1);
  dim3 dimGrid(N_E,1);
#endif

  average_cell_p0 arch_args (N_s, N_E, N_F, DU);
}

