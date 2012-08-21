#ifndef MISC_H
#define MISC_H

#include <scalar_def.h>
#include <blas_stuff.h>

void blasScopy(int N, float* x, int INCX, float* y, int INCY);
void blasSaxpy(int M, float alpha, float* x, int INCX, float* y, int INCY);
void blasSgemm(char or1, char or2, int M , int N, int K, float alpha, float* A, int LDA, float* B, int LDB, float beta, float* C, int LDC);
void blasDcopy(int N, double* x, int INCX, double* y, int INCY);
void blasDaxpy(int M, double alpha, double* x, int INCX, double* y, int INCY);
void blasDgemm(char or1, char or2, int M , int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double beta, double* C, int LDC);

void makeZero(scalar* A, int size);

#endif
