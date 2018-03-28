/*!
  \file printer_oned_singlefluid.cu
  \brief Kernel to output 1D singlefluid solution used by the PRINTER class
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup printer
*/
#ifdef ONED
#ifdef SCALARAD
#include "printer_kernels.h"
#include "constants.h"

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void formater(int N_s, int N_E, scalar* U, scalar* output, bool inverse){
  /*!
    \brief Format solution kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] U solution to format to output
    \param[out] output output solution array
    \param[in] inverse true if you want to copy from output to U instead
  */

  //
  // Copy from U to output
  //
  if(!inverse){
#ifdef USE_CPU
    for(int e = 0; e < N_E; e++){
      for(int i = 0; i < N_s; i++){
#elif USE_GPU
    int e = blockIdx.x*blkE+threadIdx.z;
    if (e < N_E){
      int i = threadIdx.x;{
#endif

	// Separate the fields
	scalar rho = U[(e*N_F+0)*N_s+i];
	
	output[(e*N_F+0)*N_s+i] = rho;

      } // loop on nodes
    }  // loop on elements
  } // if inverse

  //
  // Copy from output to U
  //
  else {
#ifdef USE_CPU
    for(int e = 0; e < N_E; e++){
      for(int i = 0; i < N_s; i++){
#elif USE_GPU
    int e = blockIdx.x*blkE+threadIdx.z;
    if (e < N_E){
      int i = threadIdx.x;{
#endif

	// Get fields from output
	scalar rho   = output[(e*N_F+0)*N_s+i];

	U[(e*N_F+0)*N_s+i] = rho;

      } // loop on nodes
    }  // loop on elements
  } // if inverse
}


//==========================================================================
//
//  Host C functions
//
//==========================================================================

extern "C"
void Lformater(int N_s, int N_E, scalar* U, scalar* output, bool inverse){
  /*!
    \brief Host C function to lauch format kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] U solution to format to output
    \param[out] output output solution array
    \param[in] inverse true if you want to copy from output to U instead (default false)
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

  formater arch_args (N_s, N_E, U, output, inverse);
};
#endif
#endif
