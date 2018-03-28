/*!
  \file BR2_tools.h
  \Header file for BR2 operations, including local lift functions and gradient calculation
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup gradient
*/
#ifndef BR2_TOOLS_H
#define BR2_TOOLS_H

#include "scalar_def.h"
#include <math.h>
#include "macros.h"
#include <stdio.h>



void Jump_Residual(int M_T, int N_s, int M_G, int N_N, int* BR2_Map, scalar* UgF, scalar* phigF, scalar* normals, scalar* Jac_trace, scalar* JF, scalar* weightF, scalar* BR2_resi, int ruk);
void Lift_Solve(int N_E, int N_N, int N_s, scalar* h_Minv, scalar* BR2_resi, scalar* BR2_poly, int ruk);
void populate_CellGrad(int N_E, int N_G, int N_s, int Chi, scalar* U, scalar* invJac,  scalar* dphig, scalar* phig, scalar* br2_poly, scalar* dUg, int ruk);
void populate_FaceGrad(int M_T, int M_G, int N_s, int N_N, int* BR2_Map, scalar* U,  scalar* inv_Jac_trace, scalar* phigF, scalar* dphigF, scalar* br2_poly, scalar* dUgF, int ruk);


#endif
