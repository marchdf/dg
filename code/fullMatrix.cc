// Gmsh - Copyright (C) 1997-2010 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#include <complex>
#include <string.h>
#include "fullMatrix.h"
#include <blas_stuff.h>
// //#if defined(_MSC_VER)
// //#define F77NAME(x) (x)
// //#endif

// #if !defined(F77NAME)
// #define F77NAME(x) (x##_)
// #endif

// // Specialisation of fullVector/Matrix operations using BLAS and LAPACK

#if defined(HAVE_BLAS)

template<> 
void fullVector<double>::axpy(const fullVector<double> &x,double alpha)
{
  int M = _r, INCX = 1, INCY = 1;
  F77NAME(daxpy)(&M, &alpha, x._data,&INCX, _data, &INCY);
}

template<> 
void fullVector<float>::axpy(const fullVector<float> &x,float alpha)
{
  int M = _r, INCX = 1, INCY = 1;
  F77NAME(saxpy)(&M, &alpha, x._data,&INCX, _data, &INCY);
}

template<> 
void fullMatrix<float>::scale(const float s)
{
  int N = _r * _c;
  int stride = 1;
  float ss = s;
  F77NAME(sscal)(&N, &ss,_data, &stride);
}

template<> 
void fullMatrix<double>::scale(const double s)
{
  int N = _r * _c;
  int stride = 1;
  double ss = s;
  F77NAME(dscal)(&N, &ss,_data, &stride);
}

template<> 
void fullMatrix<std::complex<double> >::scale(const double s)
{
  int N = _r * _c;
  int stride = 1;
  std::complex<double> ss(s, 0.);
  F77NAME(zscal)(&N, &ss,_data, &stride);
}

template<> 
void fullMatrix<double>::mult(const fullMatrix<double> &b, fullMatrix<double> &c) const
{
  int M = c.size1(), N = c.size2(), K = _c;
  int LDA = _r, LDB = b.size1(), LDC = c.size1();
  double alpha = 1., beta = 0.;
  F77NAME(dgemm)("N", "N", &M, &N, &K, &alpha, _data, &LDA, b._data, &LDB, 
                 &beta, c._data, &LDC);
}

template<> 
void fullMatrix<float>::mult(const fullMatrix<float> &b, fullMatrix<float> &c) const
{
  int M = c.size1(), N = c.size2(), K = _c;
  int LDA = _r, LDB = b.size1(), LDC = c.size1();
  float alpha = 1., beta = 0.;
  F77NAME(sgemm)("N", "N", &M, &N, &K, &alpha, _data, &LDA, b._data, &LDB, 
                 &beta, c._data, &LDC);
}

template<> 
void fullMatrix<std::complex<double> >::mult(const fullMatrix<std::complex<double> > &b, 
                                             fullMatrix<std::complex<double> > &c) const
{
  int M = c.size1(), N = c.size2(), K = _c;
  int LDA = _r, LDB = b.size1(), LDC = c.size1();
  std::complex<double> alpha = 1., beta = 0.;
  F77NAME(zgemm)("N", "N", &M, &N, &K, &alpha, _data, &LDA, b._data, &LDB, 
                 &beta, c._data, &LDC);
}

template<> 
void fullMatrix<double>::gemm(const fullMatrix<double> &a, const fullMatrix<double> &b, 
                              double alpha, double beta)
{
  int M = size1(), N = size2(), K = a.size2();
  int LDA = a.size1(), LDB = b.size1(), LDC = size1();
  F77NAME(dgemm)("N", "N", &M, &N, &K, &alpha, a._data, &LDA, b._data, &LDB, 
                 &beta, _data, &LDC);
}

template<> 
void fullMatrix<float>::gemm(const fullMatrix<float> &a, const fullMatrix<float> &b, 
                              float alpha, float beta)
{
  int M = size1(), N = size2(), K = a.size2();
  int LDA = a.size1(), LDB = b.size1(), LDC = size1();
  F77NAME(sgemm)("N", "N", &M, &N, &K, &alpha, a._data, &LDA, b._data, &LDB, 
                 &beta, _data, &LDC);
}

template<> 
void fullMatrix<std::complex<double> >::gemm(const fullMatrix<std::complex<double> > &a, 
                                             const fullMatrix<std::complex<double> > &b, 
                                             std::complex<double> alpha, 
                                             std::complex<double> beta)
{
  int M = size1(), N = size2(), K = a.size2();
  int LDA = a.size1(), LDB = b.size1(), LDC = size1();
  F77NAME(zgemm)("N", "N", &M, &N, &K, &alpha, a._data, &LDA, b._data, &LDB, 
                 &beta, _data, &LDC);
}

template<> 
void fullMatrix<double>::mult(const fullVector<double> &x, fullVector<double> &y) const
{
  int M = _r, N = _c, LDA = _r, INCX = 1, INCY = 1;
  double alpha = 1., beta = 0.;
  F77NAME(dgemv)("N", &M, &N, &alpha, _data, &LDA, x._data, &INCX,
                 &beta, y._data, &INCY);
}

template<> 
void fullMatrix<float>::mult(const fullVector<float> &x, fullVector<float> &y) const
{
  int M = _r, N = _c, LDA = _r, INCX = 1, INCY = 1;
  float alpha = 1., beta = 0.;
  F77NAME(sgemv)("N", &M, &N, &alpha, _data, &LDA, x._data, &INCX,
                 &beta, y._data, &INCY);
}

template<> 
void fullMatrix<std::complex<double> >::mult(const fullVector<std::complex<double> > &x, 
                                             fullVector<std::complex<double> > &y) const
{
  int M = _r, N = _c, LDA = _r, INCX = 1, INCY = 1;
  std::complex<double> alpha = 1., beta = 0.;
  F77NAME(zgemv)("N", &M, &N, &alpha, _data, &LDA, x._data, &INCX,
                 &beta, y._data, &INCY);
}

#endif


// FOR FLOATS
#if defined(HAVE_LAPACK)

extern "C" {
  void F77NAME(sgesv)(int *N, int *nrhs, float *A, int *lda, int *ipiv,
                      float *b, int *ldb, int *info);
  void F77NAME(sgetrf)(int *M, int *N, float *A, int *lda, int *ipiv, int *info);
  void F77NAME(sgetri)(int *M, float *A, int *lda, int *ipiv, float *work, 
                       int *lwork, int *info);
  void F77NAME(sgesvd)(const char* jobu, const char *jobvt, int *M, int *N,
                       float *A, int *lda, float *S, float* U, int *ldu,
                       float *VT, int *ldvt, float *work, int *lwork, int *info);
  void F77NAME(sgeev)(const char *jobvl, const char *jobvr, int *n, float *a,
                      int *lda, float *wr, float *wi, float *vl, int *ldvl, 
                      float *vr, int *ldvr, float *work, int *lwork, int *info);
}

static void swap(float *a, int inca, float *b, int incb, int n)
{
  float tmp;
  for (int i = 0; i < n; i++, a += inca, b += incb) {
    tmp = (*a);
    (*a) = (*b);
    (*b) = tmp;
  }
}

void eigSort(int n, float *wr, float *wi, float *VL, float *VR)
{
  // Sort the eigenvalues/vectors in ascending order according to
  // their real part. Warning: this will screw up the ordering if we
  // have complex eigenvalues.
  for (int i = 0; i < n - 1; i++){
    int k = i;
    float ek = wr[i];
    // search for something to swap
    for (int j = i + 1; j < n; j++){
      const float ej = wr[j];
      if(ej < ek){
        k = j;
        ek = ej;
      }
    }
    if (k != i){
      swap(&wr[i], 1, &wr[k], 1, 1);
      swap(&wi[i], 1, &wi[k], 1, 1);
      swap(&VL[n * i], 1, &VL[n * k], 1, n);
      swap(&VR[n * i], 1, &VR[n * k], 1, n);
    }
  }
}

template<> 
bool fullMatrix<float>::eig(fullVector<float> &DR, fullVector<float> &DI,
                             fullMatrix<float> &VL, fullMatrix<float> &VR,
                             bool sortRealPart)
{
  int N = size1(), info;
  int lwork = 10 * N;
  float *work = new float[lwork];
  F77NAME(sgeev)("V", "V", &N, _data, &N, DR._data, DI._data,
                 VL._data, &N, VR._data, &N, work, &lwork, &info);
  delete [] work;

  if(info > 0)
    printf("QR Algorithm failed to compute all the eigenvalues", info, info);
  else if(info < 0)
    printf("Wrong %d-th argument in eig", -info);
  else if(sortRealPart) 
    eigSort(N, DR._data, DI._data, VL._data, VR._data);
  
  return true;
}

template<> 
bool fullMatrix<float>::luSolve(const fullVector<float> &rhs, fullVector<float> &result)
{
  int N = size1(), nrhs = 1, lda = N, ldb = N, info;
  int *ipiv = new int[N];
  for(int i = 0; i < N; i++) result(i) = rhs(i);
  F77NAME(sgesv)(&N, &nrhs, _data, &lda, ipiv, result._data, &ldb, &info);
  delete [] ipiv;
  if(info == 0) return true;
  if(info > 0)
    printf("U(%d,%d)=0 in LU decomposition", info, info);
  else
    printf("Wrong %d-th argument in LU decomposition", -info);
  return false;
}

template<>
bool fullMatrix<float>::invert(fullMatrix<float> &result) const
{
  int M = size1(), N = size2(), lda = size1(), info;
  int *ipiv = new int[std::min(M, N)];
  result = *this;
  F77NAME(sgetrf)(&M, &N, result._data, &lda, ipiv, &info);
  if(info == 0){
    int lwork = M * 4;
    float *work = new float[lwork];
    F77NAME(sgetri)(&M, result._data, &lda, ipiv, work, &lwork, &info);
    delete [] work;
  }
  delete [] ipiv;
  if(info == 0) return true;
  else if(info > 0)
    printf("U(%d,%d)=0 in matrix inversion\n", info, info);
  else
    printf("Wrong %d-th argument in matrix inversion", -info);
  return false;
}

template<> 
bool fullMatrix<float>::invertInPlace()
{
  int N = size1(), nrhs = N, lda = N, ldb = N, info;
  int *ipiv = new int[N];
  float * invA = new float[N*N];

  for (int i = 0; i < N * N; i++) invA[i] = 0.;
  for (int i = 0; i < N; i++) invA[i * N + i] = 1.;

  F77NAME(sgesv)(&N, &nrhs, _data, &lda, ipiv, invA, &ldb, &info);
  memcpy(_data, invA, N * N * sizeof(float));

  delete [] invA;
  delete [] ipiv;

  if(info == 0) return true;
  if(info > 0)
    printf("U(%d,%d)=0 in matrix in place inversion", info, info);
  else
    printf("Wrong %d-th argument in matrix inversion", -info);
  return false;
}

template<> 
float fullMatrix<float>::determinant() const
{
  fullMatrix<float> tmp(*this);
  int M = size1(), N = size2(), lda = size1(), info;
  int *ipiv = new int[std::min(M, N)];
  F77NAME(sgetrf)(&M, &N, tmp._data, &lda, ipiv, &info);
  float det = 1.;
  if(info == 0){
    for(int i = 0; i < size1(); i++){
      det *= tmp(i, i);
      if(ipiv[i] != i + 1) det = -det;
    }
  }
  else if(info > 0)
    det = 0.;
  else
    printf("Wrong %d-th argument in matrix factorization", -info);
  delete [] ipiv;  
  return det;
}

template<> 
bool fullMatrix<float>::svd(fullMatrix<float> &V, fullVector<float> &S)
{
  fullMatrix<float> VT(V.size2(), V.size1());
  int M = size1(), N = size2(), LDA = size1(), LDVT = VT.size1(), info;
  int lwork = std::max(3 * std::min(M, N) + std::max(M, N), 5 * std::min(M, N));
  fullVector<float> WORK(lwork);
  F77NAME(sgesvd)("O", "A", &M, &N, _data, &LDA, S._data, _data, &LDA,
                  VT._data, &LDVT, WORK._data, &lwork, &info);
  V = VT.transpose();
  if(info == 0) return true;
  if(info > 0)
    printf("SVD did not converge");
  else
    printf("Wrong %d-th argument in SVD decomposition", -info);
  return false;
}

#endif


// FOR DOUBLES
#if defined(HAVE_LAPACK)

extern "C" {
  void F77NAME(dgesv)(int *N, int *nrhs, double *A, int *lda, int *ipiv,
                      double *b, int *ldb, int *info);
  void F77NAME(dgetrf)(int *M, int *N, double *A, int *lda, int *ipiv, int *info);
  void F77NAME(dgetri)(int *M, double *A, int *lda, int *ipiv, double *work, 
                       int *lwork, int *info);
  void F77NAME(dgesvd)(const char* jobu, const char *jobvt, int *M, int *N,
                       double *A, int *lda, double *S, double* U, int *ldu,
                       double *VT, int *ldvt, double *work, int *lwork, int *info);
  void F77NAME(dgeev)(const char *jobvl, const char *jobvr, int *n, double *a,
                      int *lda, double *wr, double *wi, double *vl, int *ldvl, 
                      double *vr, int *ldvr, double *work, int *lwork, int *info);
}

static void swap(double *a, int inca, double *b, int incb, int n)
{
  double tmp;
  for (int i = 0; i < n; i++, a += inca, b += incb) {
    tmp = (*a);
    (*a) = (*b);
    (*b) = tmp;
  }
}

void eigSort(int n, double *wr, double *wi, double *VL, double *VR)
{
  // Sort the eigenvalues/vectors in ascending order according to
  // their real part. Warning: this will screw up the ordering if we
  // have complex eigenvalues.
  for (int i = 0; i < n - 1; i++){
    int k = i;
    double ek = wr[i];
    // search for something to swap
    for (int j = i + 1; j < n; j++){
      const double ej = wr[j];
      if(ej < ek){
        k = j;
        ek = ej;
      }
    }
    if (k != i){
      swap(&wr[i], 1, &wr[k], 1, 1);
      swap(&wi[i], 1, &wi[k], 1, 1);
      swap(&VL[n * i], 1, &VL[n * k], 1, n);
      swap(&VR[n * i], 1, &VR[n * k], 1, n);
    }
  }
}

template<> 
bool fullMatrix<double>::eig(fullVector<double> &DR, fullVector<double> &DI,
                             fullMatrix<double> &VL, fullMatrix<double> &VR,
                             bool sortRealPart)
{
  int N = size1(), info;
  int lwork = 10 * N;
  double *work = new double[lwork];
  F77NAME(dgeev)("V", "V", &N, _data, &N, DR._data, DI._data,
                 VL._data, &N, VR._data, &N, work, &lwork, &info);
  delete [] work;

  if(info > 0)
    printf("QR Algorithm failed to compute all the eigenvalues", info, info);
  else if(info < 0)
    printf("Wrong %d-th argument in eig", -info);
  else if(sortRealPart) 
    eigSort(N, DR._data, DI._data, VL._data, VR._data);
  
  return true;
}

template<> 
bool fullMatrix<double>::luSolve(const fullVector<double> &rhs, fullVector<double> &result)
{
  int N = size1(), nrhs = 1, lda = N, ldb = N, info;
  int *ipiv = new int[N];
  for(int i = 0; i < N; i++) result(i) = rhs(i);
  F77NAME(dgesv)(&N, &nrhs, _data, &lda, ipiv, result._data, &ldb, &info);
  delete [] ipiv;
  if(info == 0) return true;
  if(info > 0)
    printf("U(%d,%d)=0 in LU decomposition", info, info);
  else
    printf("Wrong %d-th argument in LU decomposition", -info);
  return false;
}

template<>
bool fullMatrix<double>::invert(fullMatrix<double> &result) const
{
  int M = size1(), N = size2(), lda = size1(), info;
  int *ipiv = new int[std::min(M, N)];
  result = *this;
  F77NAME(dgetrf)(&M, &N, result._data, &lda, ipiv, &info);
  if(info == 0){
    int lwork = M * 4;
    double *work = new double[lwork];
    F77NAME(dgetri)(&M, result._data, &lda, ipiv, work, &lwork, &info);
    delete [] work;
  }
  delete [] ipiv;
  if(info == 0) return true;
  else if(info > 0)
    printf("U(%d,%d)=0 in matrix inversion\n", info, info);
  else
    printf("Wrong %d-th argument in matrix inversion", -info);
  return false;
}

template<> 
bool fullMatrix<double>::invertInPlace()
{
  int N = size1(), nrhs = N, lda = N, ldb = N, info;
  int *ipiv = new int[N];
  double * invA = new double[N*N];

  for (int i = 0; i < N * N; i++) invA[i] = 0.;
  for (int i = 0; i < N; i++) invA[i * N + i] = 1.;

  F77NAME(dgesv)(&N, &nrhs, _data, &lda, ipiv, invA, &ldb, &info);
  memcpy(_data, invA, N * N * sizeof(double));

  delete [] invA;
  delete [] ipiv;

  if(info == 0) return true;
  if(info > 0)
    printf("U(%d,%d)=0 in matrix in place inversion", info, info);
  else
    printf("Wrong %d-th argument in matrix inversion", -info);
  return false;
}

template<> 
double fullMatrix<double>::determinant() const
{
  fullMatrix<double> tmp(*this);
  int M = size1(), N = size2(), lda = size1(), info;
  int *ipiv = new int[std::min(M, N)];
  F77NAME(dgetrf)(&M, &N, tmp._data, &lda, ipiv, &info);
  double det = 1.;
  if(info == 0){
    for(int i = 0; i < size1(); i++){
      det *= tmp(i, i);
      if(ipiv[i] != i + 1) det = -det;
    }
  }
  else if(info > 0)
    det = 0.;
  else
    printf("Wrong %d-th argument in matrix factorization", -info);
  delete [] ipiv;  
  return det;
}

template<> 
bool fullMatrix<double>::svd(fullMatrix<double> &V, fullVector<double> &S)
{
  fullMatrix<double> VT(V.size2(), V.size1());
  int M = size1(), N = size2(), LDA = size1(), LDVT = VT.size1(), info;
  int lwork = std::max(3 * std::min(M, N) + std::max(M, N), 5 * std::min(M, N));
  fullVector<double> WORK(lwork);
  F77NAME(dgesvd)("O", "A", &M, &N, _data, &LDA, S._data, _data, &LDA,
                  VT._data, &LDVT, WORK._data, &lwork, &info);
  V = VT.transpose();
  if(info == 0) return true;
  if(info > 0)
    printf("SVD did not converge");
  else
    printf("Wrong %d-th argument in SVD decomposition", -info);
  return false;
}

#endif

