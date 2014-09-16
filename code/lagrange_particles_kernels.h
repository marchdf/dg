/*!
  \file lagrange_particles_kernels.h
  \brief Functions to launch LAGRANGE_PARTICLES kernels
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup lagrange_particles
*/
#ifndef LAGRANGE_PARTICLES_KERNELS_H
#define LAGRANGE_PARTICLES_KERNELS_H

#include <scalar_def.h>
#include <math.h>
#include "fullMatrix.h"
#include <macros.h>

int Lget_element_belong(scalar* position, int prev_el, const int* neighbors, const fullMatrix<scalar> XYZNodes, const int nvert, const int N_N, const int N_E);
void Lget_velocity_at_position(scalar* position, int el, int N_s, scalar* U, scalar* solution, scalar* avg_velocity);

#endif
