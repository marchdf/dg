#ifndef BLAS_STUFF_H
#define BLAS_STUFF_H

#include <complex>
#include <string.h>

//#if defined(_MSC_VER)
//#define F77NAME(x) (x)
//#endif

#if !defined(F77NAME)
#define F77NAME(x) (x##_)
#endif

// Specialisation of fullVector/Matrix operations using BLAS and LAPACK

#if defined(HAVE_BLAS)

extern "C" {
  void F77NAME(dcopy)(int *n, double *x, int *incx, double *y, int *incy);
  void F77NAME(scopy)(int *n, float *x, int *incx, float *y, int *incy);
  void F77NAME(daxpy)(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
  void F77NAME(saxpy)(int *n, float  *alpha, float  *x, int *incx, float  *y, int *incy);
  void F77NAME(dgemm)(const char *transa, const char *transb, int *m, int *n, int *k, 
                      double *alpha, double *a, int *lda, 
                      double *b, int *ldb, double *beta, 
                      double *c, int *ldc);
  void F77NAME(sgemm)(const char *transa, const char *transb, int *m, int *n, int *k, 
                      float *alpha, float *a, int *lda, 
                      float *b, int *ldb, float *beta, 
                      float *c, int *ldc);
  void F77NAME(zgemm)(const char *transa, const char *transb, int *m, int *n, int *k, 
                      std::complex<double> *alpha, std::complex<double> *a, int *lda, 
                      std::complex<double> *b, int *ldb, std::complex<double> *beta, 
                      std::complex<double> *c, int *ldc);
  void F77NAME(dgemv)(const char *trans, int *m, int *n, 
                      double *alpha, double *a, int *lda, 
                      double *x, int *incx, double *beta, 
                      double *y, int *incy);
  void F77NAME(sgemv)(const char *trans, int *m, int *n, 
                      float *alpha, float *a, int *lda, 
                      float *x, int *incx, float *beta, 
                      float *y, int *incy);
  void F77NAME(zgemv)(const char *trans, int *m, int *n, 
                      std::complex<double> *alpha, std::complex<double> *a, int *lda, 
                      std::complex<double> *x, int *incx, std::complex<double> *beta, 
                      std::complex<double> *y, int *incy);
  void F77NAME(dscal)(int *n, double *alpha, double *x,int *incx);
  void F77NAME(sscal)(int *n, float *alpha,  float *x, int *incx);
  void F77NAME(zscal)(int *n, std::complex<double> *alpha,std::complex<double> *x,  int *incx);
}
#endif
#endif
