/*!
  \file kernels.h
  \brief Functions to launch some kernels
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
*/
#ifndef KERNELS_H
#define KERNELS_H
#include <scalar_def.h>
#include <math.h>
#include <macros.h>

// Here I define the gpu kernels I will be using
extern "C" void LmapToFace(int M_s, int M_T,  int N_s, int* map, scalar* U, scalar* UF);
extern "C" void LmapToElement(int N_s, int N_E,  int M_s, int N_N, int* invmap, scalar* Q, scalar* q);
extern "C" void Lredistribute_sf(int N_G, int N_E,  scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac);
extern "C" void Lredistribute_q(int M_G, int M_T,  scalar* qJ, scalar* q, scalar* JF);
extern "C" void LaddSFQ(int N_s, int N_E,  scalar* A, scalar* S, scalar* F, scalar* Q); // A = S+F+Q
extern "C" void Lhrl1D(int N_s, int N_E, int N_G, int Nfields, int N_N, int slicenum, int* neighbors, int offxy, scalar* weight, scalar* V, scalar* A, scalar* Alim);
extern "C" void Lhrl2D(int N_s, int N_E, int N_G, int N_N, int order, scalar* XYZCen, scalar* powersXYZG, int* neighbors, int* TaylorDxIdx, int* TaylorDyIdx, scalar* weight, scalar refArea, scalar* A, scalar* Alim);
extern "C" void LChangeBasis(int size1, int size2, int N_E,  scalar* Transform, scalar* U, scalar* Unew);
#endif
