#ifndef MACROS_H
#define MACROS_H


//
// Architecture specific macros
//
#ifdef USE_CPU
#define arch_global
#define arch_device
#define arch_args
#define arch_args_array(x) 
#define arch(x)      h_ ## x
#define del(x)       delete[] x
#elif USE_GPU
#define arch_global __global__
#define arch_device __device__
#define arch_args   <<<dimGrid,dimBlock>>>
#define arch_args_array(x)   <<<dimGrid,dimBlock,x>>>
#define arch(x)      d_ ## x
#define del(x)       CUDA_SAFE_CALL(cudaFree(x));
#endif

// 
// Define some extra macros
//
#define MIN(X,Y)  ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y)  ((X) > (Y) ? (X) : (Y))
#define blkE 4 // number of elements per block on GPU
#define blkT 8 // number of faces per block on GPU

#endif
