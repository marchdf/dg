/*!
  \file misc_cuda.h
  \brief Small functions that run on host but that interface with the GPU.
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifdef USE_GPU
#ifndef MISC_CUDA_H
#define MISC_CUDA_H
#include <stdio.h>
#include <assert.h>
#include <cuda_runtime_api.h>

inline cudaError_t checkCuda(cudaError_t result){
  /*!
    \brief Check cuda function execution.
    \section Description
    Convenience function for checking CUDA runtime API results can be
    wrapped around any runtime API call. No-op in release builds.
    From http://devblogs.nvidia.com/parallelforall/how-optimize-data-transfers-cuda-cc/
  */
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %sn",
	    cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
  return result;
}
#endif
#endif
