/*!
  \file sensor_kernels.h
  \brief Functions to launch SENSOR kernels
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup sensor
*/
#ifndef SENSOR_KERNELS_H
#define SENSOR_KERNELS_H

#include <scalar_def.h>
#include <math.h>
#include <macros.h>

extern "C" void Lsensor1(int N_s, int N_E, scalar* U, int* sensors);
extern "C" void Lcopy_detected(int N_s, int N_E, int* sensors, scalar* Uold, scalar* U);

#endif
