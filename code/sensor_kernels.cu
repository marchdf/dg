/*!
  \file sensor_kernels.cu
  \brief Kernels used by the SENSOR class
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup sensor
*/
#include <sensor_kernels.h>


//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void sensor1(int N_s, int N_E, scalar* U, int* sensors){
  /*!
    \brief sensor1 kernel
    \param[in] N_s number of nodes per element.
    \param[in] N_E number of elements
    \param[in] U main solution
    \param[out] sensors array to hold the sensor
  */
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
#endif
    sensors[e] = 1;    
  }
}

arch_global void copy_detected(int N_s, int N_E, int* sensors, scalar* Uold, scalar* U){
  /*!
    \brief Copy only the elements which were detected with a sensor and limited.
    \param[in] N_s number of nodes per element.
    \param[in] N_E number of elements
    \param[in] sensors array to hold the sensor
    \param[in] Uold solution containing the limited solutions
    \param[out] U solution to receive the limited elements
    \section Description
    In GPU mode, launches N_E blocks of N_s x N_F x 1 threads. 
  */

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    int sen = sensors[e];
    if (sen != 0){
      for(int i = 0; i < N_s; i++){
	for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int e = blockIdx.x;
  if (e < N_E){
    int sen = sensors[e];
    if (sen != 0){
      int i = threadIdx.x;{
	int fc= threadIdx.y;{
#endif  
	  U[(e*N_F+fc)*N_s+i] = Uold[(e*N_F+fc)*N_s+i];
	} // loop on fields
      } // loop on nodes
    } // if condition on sensor
  } // loop on element
}

  
//==========================================================================
//
//  Host C functions
//
//==========================================================================
extern "C"
void Lsensor1(int N_s, int N_E, scalar* U, int* sensors){
  /*!
    \brief Host C function to launch sensor1 kernel.
    \param[in] N_s number of nodes per element.
    \param[in] N_E number of elements
    \param[in] U main solution
    \param[out] sensors array to hold the sensor
    \section Description
    In GPU mode, launches N_E/blkE blocks of 1 x 1 x blkE
    threads. blkE controls the number of elements to set on each block
  */

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(1,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

sensor1 arch_args (N_s, N_E, U, sensors);
};

extern "C"
void Lcopy_detected(int N_s, int N_E, int* sensors, scalar* Uold, scalar* U){
  /*!
    \brief Host C function to launch copy_detected kernel.
    \param[in] N_s number of nodes per element.
    \param[in] N_E number of elements
    \param[in] sensors array to hold the sensor
    \param[in] Uold solution containing the limited solutions
    \param[out] U solution to receive the limited elements
    \section Description
    In GPU mode, launches N_E blocks of N_s x N_F x 1 threads. 
  */
#ifdef USE_GPU
  dim3 dimBlock(N_s,N_F,1);
  dim3 dimGrid(N_E,1);
#endif

  copy_detected arch_args (N_s, N_E, sensors, Uold, U);
};
