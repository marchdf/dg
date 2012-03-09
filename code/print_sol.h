#ifndef PRINT_SOL_H
#define PRINT_SOL_H
#include "fullMatrix.h"
#include "simpleMesh.h"
#include <scalar_def.h>
#include <stdlib.h>

// Define the different initial condition functions
void print_dg_shallow(const int N_s, const int N_E, const int N_F, const fullMatrix<scalar> &U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append);
void print_dg_shallow(const int N_s, const int N_E, const int N_F, scalar* U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append);

void print_dg_mhd(const int N_s, const int N_E, const int N_F, const fullMatrix<scalar> &U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append, const int all, const scalar gamma);
void print_dg_mhd(const int N_s, const int N_E, const int N_F, scalar* U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append, const int all, const scalar gamma);

void print_dg_multifluid(const int N_s, const int N_E, const int N_F, const int model, const fullMatrix<scalar> &U, const simpleMesh m, const int msh_lin, const int step, const double time, const int append, const int all);
void print_dg_multifluid(const int N_s, const int N_E, const int N_F, const int model, scalar* U, const simpleMesh m, const int msh_lin, const int step, const double time, const int append, const int all);
void print_dg_multifluid_err(const int N_s, const int N_E, const int N_F, const int model, scalar* U, const simpleMesh m, const int msh_lin, const int all);

void print_dg_passive(const int N_s, const int N_E, const int N_F, scalar gamma, const fullMatrix<scalar> &U, const simpleMesh m, const int msh_lin, const int step, const double time, const int append, const int all);
void print_dg_passive(const int N_s, const int N_E, const int N_F, scalar gamma, scalar* U, const simpleMesh m, const int msh_lin, const int step, const double time, const int append, const int all);
void print_dg_passive_err(const int N_s, const int N_E, const int N_F, scalar gamma, scalar* U, const simpleMesh m, const int msh_lin, const int all);


#endif
