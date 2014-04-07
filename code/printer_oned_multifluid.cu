/*!
  \file printer_oned_multifluid.cu
  \brief Kernel to output 1D multifluid solution used by the PRINTER class
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup printer
*/
#ifdef ONED
#ifdef MULTIFLUID
#include <printer_kernels.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void formater(int N_s, int N_E, scalar* U, scalar* output){
  /*!
    \brief Format solution kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] U solution to format to output
    \param[out] output output solution array
  */
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;{
  if (e < N_E){
    int i = threadIdx.x;
#endif

    // Separate the fields
    scalar rho = U[(e*N_F+0)*N_s+i];
    scalar ux  = U[(e*N_F+1)*N_s+i]/rho;
    scalar et  = U[(e*N_F+2)*N_s+i];
#ifdef GAMCONS
    scalar gamma = 1+rho/U[(e*N_F+3)*N_s+i];
#elif  GAMNCON
    scalar gamma = 1+1.0/U[(e*N_F+3)*N_s+i];
#endif

    output[(e*N_F+0)*N_s+i] = rho;
    output[(e*N_F+1)*N_s+i] = ux;
    output[(e*N_F+2)*N_s+i] = gamma;
    output[(e*N_F+3)*N_s+i] = (gamma-1)*(et - 0.5*ux*ux*rho);
      
    // Mass fractions
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) output[(e*N_F+4+x)*N_s+i] = U[(e*N_F+4+x)*N_s+i]/rho;
#include "loop.h"
  }
  }
}


//==========================================================================
//
//  Host C functions
//
//==========================================================================

extern "C"
void Lformater(int N_s, int N_E, scalar* U, scalar* output){
  /*!
    \brief Host C function to lauch format kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] U solution to format to output
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

  formater arch_args (N_s, N_E, U, output);
};
#endif
#endif
