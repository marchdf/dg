#ifndef PRINT_SOL_H
#define PRINT_SOL_H
#include "fullMatrix.h"
#include "simpleMesh.h"
#include <scalar_def.h>
#include <stdlib.h>
#include <constants.h>
  

void print_dg(const int N_s, const int N_E, scalar* U, const simpleMesh m, const int elem_type, const int step, const double time, const int append);

#endif


/* // Define the different initial condition functions */
/* void print_dg_shallow(const int N_s, const int N_E, const fullMatrix<scalar> &U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append); */
/* void print_dg_shallow(const int N_s, const int N_E, scalar* U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append); */

/* void print_dg_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append, const int all, const scalar gamma); */
/* void print_dg_mhd(const int N_s, const int N_E, scalar* U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append, const int all, const scalar gamma); */
