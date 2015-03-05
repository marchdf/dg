/*!
  \file physics.h
  \brief Functions to lauch kernels that deal with the physics
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef PHYSICS_H
#define PHYSICS_H
#include "scalar_def.h"
#include <math.h>
#include <macros.h>

// Generic
extern "C" void Levaluate_sf(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac);//, scalar* xyz);
extern "C" void Levaluate_q(int M_G, int M_T, scalar* q, scalar* UgF, scalar* normals);//, scalar* xyzf);

extern "C" void Lkinetic_energy1D(int N_s, int N_E, scalar* rho, scalar* rhou, scalar* K);
extern "C" void Lkinetic_energy2D(int N_s, int N_E, scalar* rho, scalar* rhou, scalar* rhov, scalar* K);
extern "C" void Lpressure(int N_s, int N_E, scalar* U, scalar* p);

/* // Possibly broken: */
/* extern "C" void Levaluate_sf_shallow(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0); */
/* extern "C" void Levaluate_sf_mhd(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma); */
/* extern "C" void Levaluate_q_shallow(int M_G, int M_T, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals); */
/* extern "C" void Levaluate_q_mhd(int M_G, int M_T, scalar* q, scalar* UgF, scalar gamma, scalar* normals); */


#endif
