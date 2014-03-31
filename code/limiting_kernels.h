/*!
  \file limiting_kernels.h
  \brief Functions to launch Limiting kernels
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup limiting
*/
#ifndef LIMITING_KERNELS_H
#define LIMITING_KERNELS_H
#include <scalar_def.h>
#include <math.h>
#include <macros.h>

// Strided copy (cpu and gpu version) used in limiting procedure
extern "C" void Lstridedcopy(int count, int blocklen, int strideA, int strideB, int offsetA, int offsetB, scalar* A, scalar* B);

// Reconstruct the monomial energy coefficients with the internal and kinetic energies
extern "C" void Lreconstruct_energy(int N_s, int N_E, int slicenum, scalar* rhoeLim, scalar* KLim, scalar* EMono, scalar* ELim);

// Reconstruct the internal energy monomial coefficients using pressure and gamma in a non-oscillatory fashion
extern "C" void Linternal_energy_multifluid(int N_s, int N_E, int slicenum, scalar* p, scalar* g, scalar* rhoe);
extern "C" void Linternal_energy_stiffened(int N_s, int N_E, int slicenum, scalar* p, scalar* g, scalar* b, scalar* rhoe);

extern "C" void Lhrl1D(int N_s, int N_E, int Nfields, int N_N, int slicenum, int* neighbors, int offxy, scalar* A, scalar* Alim);
extern "C" void Lhri1D(int N_s, int N_E, int N_N, int* neighbors, int N_s1D, int slicenum, int offxy, scalar* Lag2Mono, scalar* Mono2Lag, scalar* U);
extern "C" void Lm2i1D(int N_s, int N_E, int N_N, int* neighbors, int N_s1D, int slicenum, int offxy, scalar* Lag2Mono, scalar* Mono2Lag, scalar* U);
extern "C" void Lhrl2D(int N_s, int N_E, int N_G, int N_N, int order, scalar* XYZCen, scalar* powersXYZG, int* neighbors, int* TaylorDxIdx, int* TaylorDyIdx, scalar* weight, scalar refArea, scalar* A, scalar* Alim);
extern "C" void LChangeBasis(int size1, int size2, int N_E,  scalar* Transform, scalar* U, scalar* Unew);

#endif
