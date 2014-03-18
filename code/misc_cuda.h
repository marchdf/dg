#ifdef USE_GPU
#ifndef MISC_CUDA_H
#define MISC_CUDA_H
// Small functions that run on host but that interface with the GPU.
// It requires some cuda stuff
#include <assert.h>
#include <cuda_runtime_api.h>

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
// From http://devblogs.nvidia.com/parallelforall/how-optimize-data-transfers-cuda-cc/
inline
cudaError_t checkCuda(cudaError_t result){
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %sn",
	    cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
  return result;
}
#endif
#endif
