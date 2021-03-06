/*!
  \file misc.h
  \brief Header for miscellaneous functions used by host functions
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef MISC_H
#define MISC_H

#include "scalar_def.h"
#include "blas_stuff.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "fullMatrix.h"
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include "macros.h"
#include <map>

void blasScopy(int N, float* x, int INCX, float* y, int INCY);
void blasSaxpy(int M, float alpha, float* x, int INCX, float* y, int INCY);
void blasSgemm(char or1, char or2, int M , int N, int K, float alpha, float* A, int LDA, float* B, int LDB, float beta, float* C, int LDC);
int blasIsamax(int N, float* x, int incx);
void blasDcopy(int N, double* x, int INCX, double* y, int INCY);
void blasDaxpy(int M, double alpha, double* x, int INCX, double* y, int INCY);
void blasDgemm(char or1, char or2, int M , int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double beta, double* C, int LDC);
int blasIdamax(int N, double* x, int INCX);

void makeZero(scalar* A, int size);

//template<typename T>
//void hostDeepCopyArray(int* src, int* dest, int size);

int factorial(int n);//{ return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}

void readTable(const char *fileName, fullMatrix<scalar> &XWUP, scalar &gamma, scalar &alpha, scalar &Q);
scalar interpolate(scalar x, std::vector<std::pair<scalar, scalar> > table, scalar BCL, scalar BCR);
void get_node_rgb_map(const char *fileName,std::map<int,std::vector<int> > &node_rgb_map);

#endif
