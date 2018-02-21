/*!
  \file error_quant.h  
  \brief Header file for error quantification functions
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef ERROR_QUANT_H
#define ERROR_QUANT_H
#include "fullMatrix.h"
#include "scalar_def.h"
#include "constants.h"
#include "misc.h"
#include "simpleMesh.h"
#include <vector>

void Exact_sinphil(scalar* XYZ, scalar* sol, scalar Ldomain);
void Exact_zero(scalar* XYZ, scalar* sol);
void Exact_normvtx(scalar* XYZ, scalar* sol);
void Exact_HiOWvtx(scalar* XYZ, scalar* sol);
#endif
