/*!
  \file macros.h
  \brief Macros used throughout the code
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
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
#define del(x)       cudaFree(x);
#endif

// 
// Define some extra macros
//
#define MIN(X,Y)  ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y)  ((X) > (Y) ? (X) : (Y))
#define blkE 32 // number of elements per block on GPU
#define blkT 256 // number of faces per block on GPU
#define blkComm 1 // number of communication elements per block on GPU
#endif
