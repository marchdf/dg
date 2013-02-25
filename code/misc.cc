#include <misc.h>

void blasScopy(int N, float* x, int INCX, float* y, int INCY){
  F77NAME(scopy)(&N, x, &INCX, y, &INCY);
}
void blasSaxpy(int M, float alpha, float* x, int INCX, float* y, int INCY){
  F77NAME(saxpy)(&M, &alpha, x ,&INCX, y, &INCY);
}
void blasSgemm(char or1, char or2, int M , int N, int K, float alpha, float* A, int LDA, float* B, int LDB, float beta, float* C, int LDC){
  F77NAME(sgemm)(&or1, &or2, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
void blasDcopy(int N, double* x, int INCX, double* y, int INCY){
  F77NAME(dcopy)(&N, x, &INCX, y, &INCY);
}
void blasDaxpy(int M, double alpha, double* x, int INCX, double* y, int INCY){
  F77NAME(daxpy)(&M, &alpha, x ,&INCX, y, &INCY);
}
void blasDgemm(char or1, char or2, int M , int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double beta, double* C, int LDC){
  F77NAME(dgemm)(&or1, &or2, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}

void makeZero(scalar* A, int size){
  for(int k=0; k < size; k++) A[k] = 0.0;
}

int factorial(int n){
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
