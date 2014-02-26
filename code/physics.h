#ifndef PHYSICS_H
#define PHYSICS_H
#include <scalar_def.h>
#include <math.h>
#include <macros.h>

// Generic
extern "C" void Levaluate_sf(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac);//, scalar* xyz);
extern "C" void Levaluate_q(int M_G, int M_T, scalar* q, scalar* UgF, scalar* normals);//, scalar* xyzf);

extern "C" void Lkinetic_energy1D(int N_s, int N_E, scalar* rho, scalar* rhou, scalar* K);
extern "C" void Lkinetic_energy2D(int N_s, int N_E, scalar* rho, scalar* rhou, scalar* rhov, scalar* K);
extern "C" void Lpressure(int N_s, int N_E, scalar* U, scalar* p);
extern "C" void Lpressure_u(int N_s, int N_E, scalar* U, scalar* p, scalar* u);
extern "C" void Llimmodif(int N_s, int N_E, int slicenum, scalar* A, scalar* plim, scalar* Alim);
extern "C" void Llimmodif2(int N_s, int N_E, scalar* A, scalar* plim, scalar* ulim, scalar* Alim);


// Possibly broken:
extern "C" void Lcpu_evaluate_sf_shallow(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0);
extern "C" void Lcpu_evaluate_sf_mhd(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma);
extern "C" void Lcpu_evaluate_q_shallow(int M_G, int M_T, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals);
extern "C" void Lcpu_evaluate_q_mhd(int M_G, int M_T, scalar* q, scalar* UgF, scalar gamma, scalar* normals);


#endif
