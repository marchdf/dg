/*!
  \file communicator_gpu_kernels.cu
  \brief Kernels for communicators
  \copyright Copyright (C) 2014, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \section Description
*/
#ifdef USE_MPI
#ifdef USE_GPU
#include <communicator_gpu_kernels.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
__global__ void packager(int N_s, int N_ghosts, int Nfields, int* ghostElementSend, scalar* buffer, scalar* U){
  /*!
    \brief Kernel to package data to a GPU buffer before communicating to CPU
    \param[in] N_s number of nodes per element
    \param[in] N_ghosts number of ghost elements
    \param[in] Nfields number of fields to operate on
    \param[in] ghostElementSend list of element indices to communicate
    \param[out] buffer buffer to hold data to communicate to CPU
    \param[in] U solution to communicate
  */    
  int k = blockIdx.x*blkComm + threadIdx.z;
  if (k<N_ghosts){
    int i = threadIdx.x;
    int fc= threadIdx.y;

    // local index of U to send
    int e = ghostElementSend[k];

    // Copy data from U to buffer
    buffer[(k*Nfields+fc)*N_s+i] = U[(e*Nfields+fc)*N_s+i];
  }  
}

//==========================================================================
__global__ void unpackager(int N_s, int N_ghosts, int Nfields, int* ghostElementRecv, scalar* buffer, scalar* U){
  /*!
    \brief Kernel to unpackage data from a GPU buffer after receiving from a CPU
    \param[in] N_s number of nodes per element
    \param[in] N_ghosts number of ghost elements
    \param[in] Nfields number of fields to operate on
    \param[in] ghostElementRecv list of element indices to communicate
    \param[in] buffer buffer to hold data to distribute to GPU
    \param[out] U solution to communicate
  */

  int k = blockIdx.x*blkComm + threadIdx.z;
  if (k<N_ghosts){
    int i = threadIdx.x;
    int fc= threadIdx.y;

    // local index of U to send
    int e = ghostElementRecv[k];

    // Copy data from buffer to U
    U[(e*Nfields+fc)*N_s+i] = buffer[(k*Nfields+fc)*N_s+i];
  }  
}



//==========================================================================
//
//  Host C functions
//
//==========================================================================
extern "C" 
void Lpackager(int N_s, int N_ghosts, int Nfields, int* ghostElementSend, scalar* buffer, scalar* U){
  /*!
    \brief Host C function to launch packager kernel
    \param[in] N_s number of nodes per element
    \param[in] N_ghosts number of ghost elements
    \param[in] Nfields number of fields to operate on
    \param[in] ghostElementSend list of element indices to communicate
    \param[out] buffer buffer to hold data to communicate to CPU
    \param[in] U solution to communicate
    \section Description
    In GPU mode, launches N_ghosts blocks of N_s x Nfields threads (if
    blkComm=1). blkComm can be set to increase the number of ghost
    elements per block.
  */
  int div = N_ghosts/blkComm;
  int mod = 0;
  if (N_ghosts%blkComm != 0) mod = 1;
  dim3 dimBlock(N_s,Nfields,blkComm);
  dim3 dimGrid(div+mod,1);

  packager<<<dimGrid,dimBlock>>> (N_s, N_ghosts, Nfields, ghostElementSend, buffer, U);
};

extern "C" 
void Lunpackager(int N_s, int N_ghosts, int Nfields, int* ghostElementRecv, scalar* buffer, scalar* U){
  /*!
    \brief Host C function to launch unpackager kernel
    \param[in] N_s number of nodes per element
    \param[in] N_ghosts number of ghost elements
    \param[in] Nfields number of fields to operate on
    \param[in] ghostElementRecv list of element indices to communicate
    \param[in] buffer buffer to hold data to distribute to GPU
    \param[out] U solution to communicate
    \section Description
    In GPU mode, launches N_ghosts blocks of N_s x Nfields threads (if
    blkComm=1). blkComm can be set to increase the number of ghost
    elements per block.
  */
  int div = N_ghosts/blkComm;
  int mod = 0;
  if (N_ghosts%blkComm != 0) mod = 1;
  dim3 dimBlock(N_s,Nfields,blkComm);
  dim3 dimGrid(div+mod,1);

  unpackager<<<dimGrid,dimBlock>>> (N_s, N_ghosts, Nfields, ghostElementRecv, buffer, U);
};

#endif
#endif
