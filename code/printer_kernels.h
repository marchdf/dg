/*!
  \file printer_kernels.h
  \brief Functions to launch PRINTER kernels
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup printer
*/
#ifndef PRINTER_KERNELS_H
#define PRINTER_KERNELS_H
#include <scalar_def.h>
#include <macros.h>

extern "C" void Lformater(int N_s, int N_E, scalar* U, scalar* output);
extern "C" void Lformat_sensor(int N_s, int N_E, int* sensork, scalar* output);

#endif
