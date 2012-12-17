#ifndef PHYSICS_H
#define PHYSICS_H
#include <scalar_def.h>
#include <math.h>
#include <macros.h>

// Generic
extern "C" void Levaluate_sf_1D(int D, int N_G, int N_E, int N_F, scalar gamma, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac);
extern "C" void Levaluate_q_1D(int M_G, int M_T, int N_F, scalar gamma, scalar* q, scalar* UgF, scalar* normals);

extern "C" void Levaluate_sf_2D(int D, int N_G, int N_E, int N_F, scalar gamma, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac);
extern "C" void Levaluate_q_2D(int M_G, int M_T, int N_F, scalar gamma, scalar* q, scalar* UgF, scalar* normals);


// Possibly broken:
extern "C" void Lcpu_evaluate_sf_shallow(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0);
extern "C" void Lcpu_evaluate_sf_mhd(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma);
extern "C" void Lcpu_evaluate_q_shallow(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals);
extern "C" void Lcpu_evaluate_q_mhd(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar gamma, scalar* normals);


#endif
