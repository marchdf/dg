#ifndef MISC_H
#define MISC_H

#include <scalar_def.h>
#include <blas_stuff.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "fullMatrix.h"
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort

void blasScopy(int N, float* x, int INCX, float* y, int INCY);
void blasSaxpy(int M, float alpha, float* x, int INCX, float* y, int INCY);
void blasSgemm(char or1, char or2, int M , int N, int K, float alpha, float* A, int LDA, float* B, int LDB, float beta, float* C, int LDC);
int blasIsamax(int N, float* x, int incx);
void blasDcopy(int N, double* x, int INCX, double* y, int INCY);
void blasDaxpy(int M, double alpha, double* x, int INCX, double* y, int INCY);
void blasDgemm(char or1, char or2, int M , int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double beta, double* C, int LDC);
int blasIdamax(int N, double* x, int INCX);

void makeZero(scalar* A, int size);

int factorial(int n);//{ return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}

void readTable(const char *fileName, fullMatrix<scalar> &XWUP, scalar &gamma, scalar &alpha);
scalar interpolate(scalar x, std::vector<std::pair<scalar, scalar> > table, scalar BCL, scalar BCR);

  
#endif
