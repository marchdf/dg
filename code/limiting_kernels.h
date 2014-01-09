#ifndef LIMITING_KERNELS_H
#define LIMITING_KERNELS_H
#include <scalar_def.h>
#include <macros.h>

// Strided copy (cpu and gpu version) used in limiting procedure
extern "C" void Lstridedcopy(int count, int blocklen, int strideA, int strideB, int offsetA, int offsetB, scalar* A, scalar* B);

// Reconstruct the monomial energy coefficients with the internal and kinetic energies
extern "C" void Lreconstruct_energy(int N_s, int N_E, int slicenum, scalar* rhoeLim, scalar* KLim, scalar* EMono, scalar* ELim);

// Reconstruct the internal energy monomial coefficients using pressure and gamma in a non-oscillatory fashion
extern "C" void Linternal_energy_multifluid(int N_s, int N_E, int slicenum, scalar* p, scalar* g, scalar* rhoe);
extern "C" void Linternal_energy_stiffened(int N_s, int N_E, int slicenum, scalar* p, scalar* g, scalar* b, scalar* rhoe);
#endif
