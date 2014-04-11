/*!
  \file printer_sensor.cu
  \brief Kernel to output the sensor
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup printer
*/
#include <printer_kernels.h>


//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void format_sensor(int N_s, int N_E, int* sensor, scalar* output){
  /*!
    \brief Format sensor kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] sensor sensor array
    \param[out] output output solution array
  */

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    int sen = sensor[e];
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;{
  if (e < N_E){
    int sen = sensor[e];
    int i = threadIdx.x;
#endif

    output[e*N_s+i] = sen;
  }
  }
}


//==========================================================================
//
//  Host C functions
//
//==========================================================================

extern "C"
void Lformat_sensor(int N_s, int N_E, int* sensor, scalar* output){
  /*!
    \brief Host C function to lauch format kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] sensor sensor array
    \param[out] output output solution array
    \section Description
    In GPU mode, launches N_E/blkE blocks of N_s x 1 x blkE
    threads. blkE controls the number of elements to set on each block
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  format_sensor arch_args (N_s, N_E, sensor, output);
};
