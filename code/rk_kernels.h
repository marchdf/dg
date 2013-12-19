#ifndef RK_KERNELS_H
#define RK_KERNELS_H
#include <scalar_def.h>
#include <macros.h>

extern "C" void Lsolver(int N_s, int N_E, scalar Dt, scalar* Minv, scalar* f, scalar* DU);
extern "C" void Laverage_cell_p0(const int N_s, const int N_E, scalar* DU);
extern "C" void LfindUPA(const int N_s, const int N_E, scalar* U, scalar* UPA);

#endif
