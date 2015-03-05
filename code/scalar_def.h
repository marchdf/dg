/*!
  \file scalar_def.h
  \brief Macros and typedef for scalar
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \section Description
  Allows the use of switching between floats and doubles at compile
  time.
*/
#ifndef SCALAR_DEF_H
#define SCALAR_DEF_H

#ifdef USE_FLOAT
typedef float scalar;
#define MPI_SCALAR MPI_FLOAT
#define hostblasGemm blasSgemm // To make sure that there is a gemm blas only on the cpu
#ifdef USE_CPU
#define blasGemm blasSgemm
#define blasAxpy blasSaxpy
#define blasCopy blasScopy
#define blasScal blasSScal
#define blasIamax blasIsamax
#elif USE_GPU
#define blasGemm cublasSgemm
#define blasAxpy cublasSaxpy
#define blasCopy cublasScopy
#define blasScal cublasSscal
#define blasIamax cublasIsamax
#endif
#endif

#ifdef USE_DOUBLE
typedef double scalar;
#define MPI_SCALAR MPI_DOUBLE
#define hostblasGemm blasDgemm  // To make sure that there is a gemm blas only on the cpu
#ifdef USE_CPU
#define blasGemm blasDgemm
#define blasAxpy blasDaxpy
#define blasCopy blasDcopy
#define blasScal blasDscal
#define blasIamax blasIdamax
#elif USE_GPU
#define blasGemm cublasDgemm
#define blasAxpy cublasDaxpy
#define blasCopy cublasDcopy
#define blasScal cublasDscal
#define blasIamax cublasIdamax
#endif
#endif

#endif
