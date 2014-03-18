#ifdef USE_MPI
#ifdef USE_GPU
#include <communicator_gpu_kernels.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
__global__ void packager(int N_s, int N_ghosts, int N_fields, int* ghostElementSend, scalar* buffer, scalar* U){
  /* Function to package data from U into a buffer so that it can be
     sent to a CPU */
    
  int k = blockIdx.x*blkComm + threadIdx.z;
  if (k<N_ghosts){
    int i = threadIdx.x;
    int fc= threadIdx.y;

    // local index of U to send
    int e = ghostElementSend[k];

    // Copy data from U to buffer
    buffer[(k*N_fields+fc)*N_s+i] = U[(e*N_fields+fc)*N_s+i];
  }  
}

//==========================================================================
__global__ void unpackager(int N_s, int N_ghosts, int N_fields, int* ghostElementRecv, scalar* buffer, scalar* U){
  /* Function to unpackage data from a buffer into U so that it can be
     received to a CPU */

  int k = blockIdx.x*blkComm + threadIdx.z;
  if (k<N_ghosts){
    int i = threadIdx.x;
    int fc= threadIdx.y;

    // local index of U to send
    int e = ghostElementRecv[k];

    // Copy data from buffer to U
    U[(e*N_fields+fc)*N_s+i] = buffer[(k*N_fields+fc)*N_s+i];
  }  
}



//==========================================================================
//
//  Host C functions
//
//==========================================================================
extern "C" 
void Lpackager(int N_s, int N_ghosts, int N_fields, int* ghostElementSend, scalar* buffer, scalar* U){
  int div = N_ghosts/blkComm;
  int mod = 0;
  if (N_ghosts%blkComm != 0) mod = 1;
  dim3 dimBlock(N_s,N_fields,blkComm);
  dim3 dimGrid(div+mod,1);

  packager<<<dimGrid,dimBlock>>> (N_s, N_ghosts, N_fields, ghostElementSend, buffer, U);
};

extern "C" 
void Lunpackager(int N_s, int N_ghosts, int N_fields, int* ghostElementRecv, scalar* buffer, scalar* U){
  int div = N_ghosts/blkComm;
  int mod = 0;
  if (N_ghosts%blkComm != 0) mod = 1;
  dim3 dimBlock(N_s,N_fields,blkComm);
  dim3 dimGrid(div+mod,1);

  unpackager<<<dimGrid,dimBlock>>> (N_s, N_ghosts, N_fields, ghostElementRecv, buffer, U);
};

#endif
#endif
