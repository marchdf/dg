/*!
  \file boundaries.h
  \brief Functions to launch BC kernels
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
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
