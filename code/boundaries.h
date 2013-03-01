#ifndef BOUNDARIES_KERNELS_H
#define BOUNDARIES_KERNELS_H
#include <scalar_def.h>
#include <math.h>
#include <macros.h>


extern "C" void LperiodicBoundary(int M_s, int N_F, int M_B, int* boundaryMap, int start, scalar* UF);
extern "C" void LfarfieldBoundary(int M_s, int N_F, int M_B, int* boundaryMap, int start, scalar* UF);


#endif
