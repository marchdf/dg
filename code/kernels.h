/*!
  \file kernels.h
  \brief Functions to launch some kernels
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef KERNELS_H
#define KERNELS_H
#include "scalar_def.h"
#include "macros.h"

// Here I define the gpu kernels I will be using
extern "C" void LmapToFace(int M_s, int M_T,  int N_s, int* map, scalar* U, scalar* UF);
extern "C" void LmapToElement(int N_s, int N_E,  int M_s, int N_N, int* invmap, scalar* Q, scalar* q);
extern "C" void Lredistribute_sf(int N_G, int N_E,  scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac);
extern "C" void Lredistribute_q(int M_G, int M_T,  scalar* qJ, scalar* q, scalar* JF);
extern "C" void LaddSFQ(int N_s, int N_E,  scalar* A, scalar* S, scalar* F, scalar* Q); // A = S+F+Q
#endif
