#include <boundaries.h>
#include <cstdlib>
#include <stdio.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void rflctive(int M_s, int M_B, int* boundaryMap, int start, scalar* UF){
#ifdef USE_CPU
  for(int k = 0; k < M_B; k++){
    int t = boundaryMap[start+k];
    for(int j = 0; j < M_s; j++){
#elif USE_GPU
      int t = boundaryMap[start+blockIdx.x];
      int j = threadIdx.x;
#endif

#ifdef ONED
      // velocity flip
      UF[((t*N_F+1)*2+1)*M_s+j] = -UF[((t*N_F+1)*2+1)*M_s+j];
#elif TWOD
      // velocities
      UF[((t*N_F+1)*2+1)*M_s+j] = -UF[((t*N_F+1)*2+1)*M_s+j];
      UF[((t*N_F+2)*2+1)*M_s+j] = -UF[((t*N_F+2)*2+1)*M_s+j];
#endif //on dimensions
      
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
void LrflctiveBoundary(int M_s, int M_B, int* boundaryMap, int start, scalar* UF){

#ifdef USE_GPU
  dim3 dimBlock(M_s,1,1);
  dim3 dimGrid(M_B,1);
#endif

  rflctive arch_args (M_s, M_B, boundaryMap, start, UF);
}


