/*!
  \file vis_physics.h
  \brief Functions to lauch kernels that deal with the viscous physics
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef VIS_PHYSICS_H
#define VIS_PHYSICS_H
#include "scalar_def.h"
#include "constants.h"
#include <math.h>
#include <macros.h>

// Generic
extern "C" void Levaluate_sf_vis(int N_G, int N_E, int PROC, scalar* f, scalar* Ug, scalar* dUg);
extern "C" void Levaluate_q_vis(int M_G, int M_T, int PROC, scalar* q, scalar* UhCommon, scalar* gradCommon, scalar* normals);//, scalar* xyzf);

#endif
