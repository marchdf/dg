/*!
  \file sensor_kernels.h
  \brief Functions to launch SENSOR kernels
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup sensor
*/
#ifndef SENSOR_KERNELS_H
#define SENSOR_KERNELS_H

#include <scalar_def.h>
#include <math.h>
#include <macros.h>
#include <constants.h>

extern "C" void Lcalc_sensors(int N_E, int N_N, bool sensor1, scalar thresh1, bool sensor2, scalar thresh2, bool sensor3, scalar thresh3, int* neighbors, scalar* Uavg, int* sensors);
extern "C" void Lcopy_detected(int N_s, int N_E, int* sensors, scalar* Uold, scalar* U);

#endif
