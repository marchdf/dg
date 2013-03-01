#ifndef CPU_KERNELS_H
#define CPU_KERNELS_H
#include <scalar_def.h>
#include <math.h>
#include <macros.h>

// Here I define the gpu kernels I will be using
extern "C" void Lcpu_equal(int N_s, int N_E, int N_F, scalar* A, scalar* B);
extern "C" void Lcpu_add(int N_s, int N_E, int N_F, scalar* A, scalar* B, scalar c); // A = A + c*B
extern "C" void Lcpu_mapToFace_shallow(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF);
extern "C" void Lcpu_mapToFace_mhd(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF);
extern "C" void Lcpu_mapToFace(int M_s, int M_T, int N_F, int N_s, int* map, scalar* U, scalar* UF);
extern "C" void Lcpu_mapToElement(int N_s, int N_E, int N_F, int* invmapx, scalar* Q, scalar* q);
extern "C" void Lcpu_collocationU(int D, int N_G, int N_s, int N_E, int N_F, scalar* Ug, scalar* dUg, scalar* phi, scalar* dphi, scalar* U);
extern "C" void Lcpu_collocationUF(int M_G, int M_s, int M_T, int N_F, scalar* UgF, scalar* psi, scalar* UF);
extern "C" void Lcpu_redistribute_sf(int D, int N_G, int N_E, int N_F, scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac);
extern "C" void Lcpu_gemm_sf(int D, int N_G, int N_s, int N_E, int N_F, scalar* S, scalar* F, scalar* sJ, scalar* fJ, scalar* phi_w, scalar* dphi_w);
extern "C" void Lcpu_redistribute_q(int M_G, int M_T, int N_F, scalar* qJ, scalar* q, scalar* JF);
extern "C" void Lcpu_gemm_q(int M_G, int M_s, int M_T, int N_F, scalar* Qtcj, scalar* qJ, scalar* psi_w);
extern "C" void Lcpu_addSFQ(int N_s, int N_E, int N_F, scalar* A, scalar* S, scalar* F, scalar* Q); // A = S+F+Q
extern "C" void Lcpu_hsl(int N_s, int N_E, int N_F, int boundaryMap, scalar* U, scalar* UNew);
extern "C" void Lcpu_hrl1D(int N_s, int N_E, int N_F, int N_G, int N_N, int slicenum, int* neighbors, int offxy, scalar* weight, scalar* V, scalar* A, scalar* Alim);
extern "C" void Lcpu_hrl2D(int N_s, int N_E, int N_F, int N_G, int N_N, int D, int order, scalar* XYZCen, scalar* powersXYZG, int* neighbors, int* TaylorDxIdx, int* TaylorDyIdx, scalar* weight, scalar refArea, scalar* A, scalar* Alim);
extern "C" void LChangeBasis(int size1, int size2, int N_E, int N_F, scalar* Transform, scalar* U, scalar* Unew);
extern "C" void Lcpu_Prim2Cons(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0);
extern "C" void Lcpu_Cons2Prim(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0);
extern "C" void Lcpu_Half2Cons(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0);
extern "C" void Lcpu_Cons2Half(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0);

#endif
