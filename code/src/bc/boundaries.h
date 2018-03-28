/*!
  \file boundaries.h
  \brief Functions to launch BC kernels
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup boundaries
*/
#ifndef BOUNDARIES_KERNELS_H
#define BOUNDARIES_KERNELS_H
#include "scalar_def.h"
#include <math.h>
#include "macros.h"
#include "constants.h"

extern "C" void LrflctiveBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF);
extern "C" void LnoslipBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF);
extern "C" void LnogradBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* dUgF);
extern "C" void LAnflwBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF);
extern "C" void LKJetBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF);
extern "C" void LSubOutBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF);
extern "C" void LHomoBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF);

scalar getMew(scalar mew_ref, scalar T_ref, scalar Cvis, scalar T);

#endif
