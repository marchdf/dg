#ifndef CPU_KERNELS_H
#define CPU_KERNELS_H
#include <scalar_def.h>
#include <math.h>
#include <macros.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

// Here I define the gpu kernels I will be using
extern "C" void Lcpu_mapToFace_shallow(int M_s, int M_T,  int* map, scalar* U, scalar* UF);
extern "C" void Lcpu_mapToFace_mhd(int M_s, int M_T,  int* map, scalar* U, scalar* UF);
extern "C" void Lcpu_mapToFace(int M_s, int M_T,  int N_s, int* map, scalar* U, scalar* UF);
extern "C" void Lcpu_mapToElement(int N_s, int N_E,  int M_s, int N_N, int* invmapx, scalar* Q, scalar* q);
extern "C" void Lcpu_redistribute_sf(int N_G, int N_E,  scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac);
extern "C" void Lcpu_redistribute_q(int M_G, int M_T,  scalar* qJ, scalar* q, scalar* JF);
extern "C" void Lcpu_addSFQ(int N_s, int N_E,  scalar* A, scalar* S, scalar* F, scalar* Q); // A = S+F+Q
extern "C" void Lcpu_hsl(int N_s, int N_E,  int boundaryMap, scalar* U, scalar* UNew);
extern "C" void Lcpu_hrl1D(int N_s, int N_E, int N_G, int Nfields, int N_N, int slicenum, int* neighbors, int offxy, scalar* weight, scalar* V, scalar* A, scalar* Alim);
extern "C" void Lcpu_hrl2D(int N_s, int N_E, int N_G, int N_N, int order, scalar* XYZCen, scalar* powersXYZG, int* neighbors, int* TaylorDxIdx, int* TaylorDyIdx, scalar* weight, scalar refArea, scalar* A, scalar* Alim);
extern "C" void LChangeBasis(int size1, int size2, int N_E,  scalar* Transform, scalar* U, scalar* Unew);
extern "C" void Lcpu_Prim2Cons(int N_s, int N_E,  scalar* U, bool multifluid, bool passive, int model, scalar gamma0);
extern "C" void Lcpu_Cons2Prim(int N_s, int N_E,  scalar* U, bool multifluid, bool passive, int model, scalar gamma0);
extern "C" void Lcpu_Half2Cons(int N_s, int N_E,  scalar* U, bool multifluid, bool passive, int model, scalar gamma0);
extern "C" void Lcpu_Cons2Half(int N_s, int N_E,  scalar* U, bool multifluid, bool passive, int model, scalar gamma0);


#endif
