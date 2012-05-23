#ifndef SCALAR_DEF_H
#define SCALAR_DEF_H

#ifdef USE_FLOAT
typedef float scalar;
#define hostblasGemm blasSgemm // To make sure that there is a gemm blas only on the cpu
#ifdef USE_CPU
#define blasGemm blasSgemm
#define blasAxpy blasSaxpy
#define blasCopy blasScopy
#define blasScal blasSScal
#elif USE_GPU
#define blasGemm cublasSgemm
#define blasAxpy cublasSaxpy
#define blasCopy cublasScopy
#define blasScal cublasSscal
#endif
#endif

#ifdef USE_DOUBLE
typedef double scalar;
#define hostblasGemm blasDgemm  // To make sure that there is a gemm blas only on the cpu
#ifdef USE_CPU
#define blasGemm blasDgemm
#define blasAxpy blasDaxpy
#define blasCopy blasDcopy
#define blasScal blasDscal
#elif USE_GPU
#define blasGemm cublasDgemm
#define blasAxpy cublasDaxpy
#define blasCopy cublasDcopy
#define blasScal cublasDscal
#endif
#endif

#endif
