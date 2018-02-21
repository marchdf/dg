/*!
  \file printer_kernels.h
  \brief Functions to launch PRINTER kernels
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup printer
*/
#ifndef PRINTER_KERNELS_H
#define PRINTER_KERNELS_H
#include "scalar_def.h"
#include "macros.h"
#include "stdio.h" //PEJ 10/12/2017

extern "C" void Lformater(int N_s, int N_E, scalar* U, scalar* output, bool inverse = false);
extern "C" void Lformat_sensor(int N_s, int N_E, int* sensor, scalar* output);

#endif
