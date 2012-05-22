#ifndef SCALAR_DEF_H
#define SCALAR_DEF_H

#ifdef USE_FLOAT
typedef float scalar;
#define blasGemm   blasSgemm
#define blasAxpy   blasSaxpy
#define blasCopy   blasScopy
#define cublasGemm cublasSgemm
#define cublasAxpy cublasSaxpy
#define cublasCopy cublasScopy
#define cublasScal cublasSscal
#endif

#ifdef USE_DOUBLE
typedef double scalar;
#define blasGemm   blasDgemm
#define blasAxpy   blasDaxpy
#define blasCopy   blasDcopy
#define cublasGemm cublasDgemm
#define cublasAxpy cublasDaxpy
#define cublasCopy cublasDcopy
#define cublasScal cublasDscal
#endif

#endif
