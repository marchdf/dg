#include <boundaries.h>
//#include <multifluid_rflctive.h>
#include <cstdlib>
#include <stdio.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================


//==========================================================================
arch_global void periodic(int M_s, int N_F, int M_B, int* boundaryMap, int start, scalar* UF){
#ifdef USE_CPU
  for(int t = 0; t < M_B; t++){
    int t1 = boundaryMap[start+t*2+0];
    int t2 = boundaryMap[start+t*2+1];
    for(int j = 0; j < M_s; j++){
      for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
	int t1 = boundaryMap[start+blockIdx.x*2+0];
	int t2 = boundaryMap[start+blockIdx.x*2+1];
	int j  = threadIdx.x;
	int fc = threadIdx.y;
#endif


#ifdef ONED
	UF[((t1*N_F+fc)*2+1)*M_s+j] = UF[((t2*N_F+fc)*2+0)*M_s+j]; // original buggy one (works for 1D)
#elif TWOD
	// This is fiddling bc the boundaries get rotated... 
	if      (j==0) UF[((t1*N_F+fc)*2+1)*M_s+0] = UF[((t2*N_F+fc)*2+0)*M_s+1];
	else if (j==1) UF[((t1*N_F+fc)*2+1)*M_s+1] = UF[((t2*N_F+fc)*2+0)*M_s+0];
	else           UF[((t1*N_F+fc)*2+1)*M_s+j] = UF[((t2*N_F+fc)*2+0)*M_s+M_s-1-j];
#endif

#ifdef USE_CPU
      }
    }
  }
#endif
}

//==========================================================================
arch_global void farfield(int M_s, int N_F, int M_B, int* boundaryMap, int start, scalar* UF){
#ifdef USE_CPU
  for(int t = 0; t < M_B; t++){
    int t1 = boundaryMap[start+t*2+0];
    int t2 = boundaryMap[start+t*2+1];
    for(int j = 0; j < M_s; j++){
      for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
	int t1 = boundaryMap[start+blockIdx.x*2+0];
	int t2 = boundaryMap[start+blockIdx.x*2+1];
	int j  = threadIdx.x;
	int fc = threadIdx.y;
#endif
	// No need to fiddle for farfield bc not rotating the boundaries
	UF[((t1*N_F+fc)*2+1)*M_s+j] = UF[((t2*N_F+fc)*2+0)*M_s+j];

#ifdef USE_CPU
      }
    }
  }
#endif
}

//==========================================================================
arch_global void rflctive(int M_s, int N_F, int M_B, int* boundaryMap, int start, scalar* UF){
#ifdef USE_CPU
  for(int t = 0; t < M_B; t++){
    int t1 = boundaryMap[start+t*2+0];
    int t2 = boundaryMap[start+t*2+1];
    for(int j = 0; j < M_s; j++){
#elif USE_GPU
	int t1 = boundaryMap[start+blockIdx.x*2+0];
	int t2 = boundaryMap[start+blockIdx.x*2+1];
	int j  = threadIdx.x;
	int fc = threadIdx.y;
#endif

#ifdef PASSIVE
	#ifdef ONED
	// density
	UF[((t1*N_F+0)*2+1)*M_s+j] = UF[((t2*N_F+0)*2+0)*M_s+j];
	// velocities
	UF[((t1*N_F+1)*2+1)*M_s+j] =-UF[((t2*N_F+1)*2+0)*M_s+j];
	// energy
	UF[((t1*N_F+2)*2+1)*M_s+j] = UF[((t2*N_F+2)*2+0)*M_s+j];
	// passive scalar
	UF[((t1*N_F+3)*2+1)*M_s+j] = UF[((t2*N_F+3)*2+0)*M_s+j];
	UF[((t1*N_F+4)*2+1)*M_s+j] = UF[((t2*N_F+4)*2+0)*M_s+j];
#elif TWOD

	// density
	UF[((t1*N_F+0)*2+1)*M_s+j] = UF[((t2*N_F+0)*2+0)*M_s+j];
	// velocities
	UF[((t1*N_F+1)*2+1)*M_s+j] =-UF[((t2*N_F+1)*2+0)*M_s+j];
	UF[((t1*N_F+2)*2+1)*M_s+j] =-UF[((t2*N_F+2)*2+0)*M_s+j];
	// energy
	UF[((t1*N_F+3)*2+1)*M_s+j] = UF[((t2*N_F+3)*2+0)*M_s+j];
	// passive scalar
	UF[((t1*N_F+4)*2+1)*M_s+j] = UF[((t2*N_F+4)*2+0)*M_s+j];
	UF[((t1*N_F+5)*2+1)*M_s+j] = UF[((t2*N_F+5)*2+0)*M_s+j];

#endif //on dimensions

#elif MULTIFLUID
#ifdef ONED
	//oned_multifluid_rflctive();
	// density
	UF[((t1*N_F+0)*2+1)*M_s+j] = UF[((t2*N_F+0)*2+0)*M_s+j];
	// velocities
	UF[((t1*N_F+1)*2+1)*M_s+j] =-UF[((t2*N_F+1)*2+0)*M_s+j];
	// energy
	UF[((t1*N_F+2)*2+1)*M_s+j] = UF[((t2*N_F+2)*2+0)*M_s+j];
	// gamma
	UF[((t1*N_F+3)*2+1)*M_s+j] = UF[((t2*N_F+3)*2+0)*M_s+j];
#elif TWOD
	//twod_multifluid_rflctive();

	// density
	UF[((t1*N_F+0)*2+1)*M_s+j] = UF[((t2*N_F+0)*2+0)*M_s+j];
	// velocities
	UF[((t1*N_F+1)*2+1)*M_s+j] = -UF[((t2*N_F+1)*2+0)*M_s+j];
	UF[((t1*N_F+2)*2+1)*M_s+j] = -UF[((t2*N_F+2)*2+0)*M_s+j];
	// energy
	UF[((t1*N_F+3)*2+1)*M_s+j] = UF[((t2*N_F+3)*2+0)*M_s+j];
	// gamma
	UF[((t1*N_F+4)*2+1)*M_s+j] = UF[((t2*N_F+4)*2+0)*M_s+j];

#endif //on dimensions
#endif //on model type

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
void LperiodicBoundary(int M_s, int N_F, int M_B, int* boundaryMap, int start, scalar* UF){

#ifdef USE_GPU
  dim3 dimBlock(M_s,N_F,1);
  dim3 dimGrid(M_B,1);
#endif

  periodic arch_args (M_s, N_F, M_B, boundaryMap, start, UF);
}

extern "C"
void LfarfieldBoundary(int M_s, int N_F, int M_B, int* boundaryMap, int start, scalar* UF){

#ifdef USE_GPU
  dim3 dimBlock(M_s,N_F,1);
  dim3 dimGrid(M_B,1);
#endif

  farfield arch_args (M_s, N_F, M_B, boundaryMap, start, UF);
}

extern "C"
void LrflctiveBoundary(int M_s, int N_F, int M_B, int* boundaryMap, int start, scalar* UF){

#ifdef USE_GPU
  dim3 dimBlock(M_s,N_F,1);
  dim3 dimGrid(M_B,1);
#endif

  rflctive arch_args (M_s, N_F, M_B, boundaryMap, start, UF);
}


