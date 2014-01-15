#ifndef BOUNDARIES_KERNELS_H
#define BOUNDARIES_KERNELS_H
#include <scalar_def.h>
#include <math.h>
#include <macros.h>

extern "C" void LrflctiveBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF);


#endif
