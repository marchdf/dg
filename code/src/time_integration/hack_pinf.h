/*!
  \file hack_pinf.h
  \brief Functions to launch hack pinf kernels
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup rk
*/
#ifndef HACK_PINF_H
#define HACK_PINF_H
#include "scalar_def.h"
#include "macros.h"
#include <math.h>
#include "constants.h"

extern "C" void Lhack_pinf(int N_s, int N_E, scalar* U);

#endif
