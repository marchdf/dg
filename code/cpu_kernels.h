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
extern "C" void Lcpu_mapToFace_multifluid(int M_s, int M_T, int N_F, int N_s, int boundaryMap, scalar* U, scalar* UF);
extern "C" void Lcpu_mapToFace_passive(int M_s, int M_T, int N_F, int N_s, int boundaryMap, scalar* U, scalar* UF);
extern "C" void Lcpu_boundary(int M_s, int N_F, int M_B, int* boundaryMap, scalar* UF);
extern "C" void Lcpu_mapToElement(int N_s, int N_E, int N_F, scalar* Q, scalar* q);
extern "C" void Lcpu_collocationU(int D, int N_G, int N_s, int N_E, int N_F, scalar* Ug, scalar* dUg, scalar* phi, scalar* dphi, scalar* U);
extern "C" void Lcpu_collocationUF(int M_G, int M_s, int M_T, int N_F, scalar* UgF, scalar* psi, scalar* UF);
extern "C" void Lcpu_evaluate_sf_shallow(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0);
extern "C" void Lcpu_evaluate_sf_mhd(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma);
extern "C" void Lcpu_evaluate_sf_multifluid(int D, int N_G, int N_E, int N_F, int model, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac);
extern "C" void Lcpu_evaluate_sf_passive(int D, int N_G, int N_E, int N_F, scalar gamma, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac);
extern "C" void Lcpu_evaluate_q_shallow(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals);
extern "C" void Lcpu_evaluate_q_mhd(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar gamma, scalar* normals);
extern "C" void Lcpu_evaluate_q_multifluid(int M_G, int M_T, int N_F, int flux, int model, scalar* q, scalar* UgF);
extern "C" void Lcpu_evaluate_q_passive(int M_G, int M_T, int N_F, int flux, scalar gamma, scalar* q, scalar* UgF);
extern "C" void Lcpu_redistribute_sf(int D, int N_G, int N_E, int N_F, scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac);
extern "C" void Lcpu_gemm_sf(int D, int N_G, int N_s, int N_E, int N_F, scalar* S, scalar* F, scalar* sJ, scalar* fJ, scalar* phi_w, scalar* dphi_w);
extern "C" void Lcpu_redistribute_q(int M_G, int M_T, int N_F, scalar* qJ, scalar* q);
extern "C" void Lcpu_gemm_q(int M_G, int M_s, int M_T, int N_F, scalar* Qtcj, scalar* qJ, scalar* psi_w);
extern "C" void Lcpu_solve(int N_s, int N_E, int N_F, scalar* DU, scalar* S, scalar* F, scalar* Q, scalar* Minv, scalar Dt);
extern "C" void Lcpu_addSFQ(int N_s, int N_E, int N_F, scalar* A, scalar* S, scalar* F, scalar* Q); // A = S+F+Q
extern "C" void Lcpu_average_cell_p0(const int N_s, const int N_E, const int N_F, scalar* DU);
extern "C" void Lcpu_hsl(int N_s, int N_E, int N_F, int boundaryMap, scalar* U, scalar* UNew);
extern "C" void Lcpu_hrl(int N_s, int N_E, int N_F, int N_G, int boundaryMap, scalar* weight, scalar* V, scalar* A, scalar* Alim);
extern "C" void Lcpu_Prim2Cons(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0);
extern "C" void Lcpu_Cons2Prim(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0);
#endif
