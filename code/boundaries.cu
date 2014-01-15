#include <boundaries.h>
#include <cstdlib>
#include <stdio.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void rflctive(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
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
      // normal
      scalar nx = normals[t*D+0]; 
      scalar ny = normals[t*D+1]; 
      scalar nx2 = nx*nx;
      scalar ny2 = ny*ny;
      scalar m2nxny = -2*nx*ny;
      scalar invnx2ny2 = 1.0/(nx2+ny2);
      
      // // Normal and tangential velocities
      scalar vxL = UF[((t*N_F+1)*2+0)*M_s+j];
      scalar vyL = UF[((t*N_F+2)*2+0)*M_s+j];
      scalar vxR = invnx2ny2 * ((-nx2+ny2)*vxL +    m2nxny*vyL);
      scalar vyR = invnx2ny2 * (   m2nxny *vxL + (nx2-ny2)*vyL);
      
      // velocities 
      UF[((t*N_F+1)*2+1)*M_s+j] = vxR;
      UF[((t*N_F+2)*2+1)*M_s+j] = vyR;

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
void LrflctiveBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){

#ifdef USE_GPU
  dim3 dimBlock(M_s,1,1);
  dim3 dimGrid(M_B,1);
#endif

  rflctive arch_args (M_s, M_B, boundaryMap, normals, start, UF);
}


