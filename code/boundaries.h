/*!
  \file boundaries.h
  \brief Functions to launch BC kernels
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \defgroup boundaries Boundary Conditions
  \ingroup boundaries
*/
#ifndef BOUNDARIES_KERNELS_H
#define BOUNDARIES_KERNELS_H
#include <scalar_def.h>
#include <math.h>
#include <macros.h>

extern "C" void LrflctiveBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF);


#endif
