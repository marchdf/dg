/*!
  \file printer_oned_singlefluid.cu
  \brief Kernel to output 1D singlefluid solution used by the PRINTER class
  \copyright Copyright (C) 2014, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup printer
*/
#ifdef ONED
#ifdef SINGLEFLUID
#include <printer_kernels.h>
#include <constants.h>

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
	scalar ux  = U[(e*N_F+1)*N_s+i]/rho;
	scalar et  = U[(e*N_F+2)*N_s+i];
	scalar gamma = constants::GLOBAL_GAMMA;
	
	output[(e*N_F+0)*N_s+i] = rho;
	output[(e*N_F+1)*N_s+i] = ux;
	output[(e*N_F+2)*N_s+i] = (gamma-1)*(et - 0.5*ux*ux*rho);

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
	scalar ux    = output[(e*N_F+1)*N_s+i];
	scalar gamma = constants::GLOBAL_GAMMA;
	scalar p     = output[(e*N_F+2)*N_s+i];

	U[(e*N_F+0)*N_s+i] = rho;
	U[(e*N_F+1)*N_s+i] = rho*ux;
	U[(e*N_F+2)*N_s+i] = p/(gamma-1) + 0.5*rho*ux*ux;

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
