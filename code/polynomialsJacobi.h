#ifndef POLYNOMIALS_JACOBI_H
#define POLYNOMIALS_JACOBI_H
#include "fullMatrix.h"
#include <scalar_def.h>
#include <stdlib.h>
#include <math.h>


void JacobiP(const fullMatrix<scalar> x, const int alpha, const int beta, const int N, fullMatrix<scalar> &P);
void JacobiGQ(const int alpha, const int beta, const int N, fullMatrix<scalar> &x, fullMatrix<scalar> &w);
void JacobiGL(const int alpha, const int beta, const int N, fullMatrix<scalar> &x);

#endif
