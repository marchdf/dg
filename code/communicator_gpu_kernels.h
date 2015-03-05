/*!
  \file communicator_gpu_kernels.h
  \brief Functions to launch communicator kernels
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifdef USE_MPI
#ifdef USE_GPU
#ifndef COMMUNICATOR_GPU_KERNELS_H
#define COMMUNICATOR_GPU_KERNELS_H
#include <scalar_def.h>
#include <macros.h>

extern "C" void   Lpackager(int N_s, int N_ghosts, int Nfields, int* ghostElementSend, scalar* buffer, scalar* U);
extern "C" void Lunpackager(int N_s, int N_ghosts, int Nfields, int* ghostElementRecv, scalar* buffer, scalar* U);

#endif
#endif
#endif
