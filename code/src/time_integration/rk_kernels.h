/*!
  \file rk_kernels.h
  \brief Functions to launch RK kernels
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup rk
*/
#ifndef RK_KERNELS_H
#define RK_KERNELS_H
#include "scalar_def.h"
#include "macros.h"

extern "C" void Lsolver(int N_s, int N_E, scalar Dt, scalar* Minv, scalar* f, scalar* DU);
extern "C" void Laverage_cell_p0(const int N_s, const int N_E, scalar* DU);
extern "C" void LfindUPA(const int N_s, const int N_E, scalar* U, scalar* UPA);

#endif
