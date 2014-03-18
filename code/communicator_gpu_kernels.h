#ifdef USE_MPI
#ifdef USE_GPU
#ifndef COMMUNICATOR_GPU_KERNELS_H
#define COMMUNICATOR_GPU_KERNELS_H
#include <scalar_def.h>
#include <macros.h>

// package data to a GPU buffer before communicating to CPU
extern "C" void   Lpackager(int N_s, int N_ghosts, int N_fields, int* ghostElementSend, scalar* buffer, scalar* U);

// unpackage data from a GPU buffer after receiving from a CPU
extern "C" void Lunpackager(int N_s, int N_ghosts, int N_fields, int* ghostElementRecv, scalar* buffer, scalar* U);

#endif
#endif
#endif
