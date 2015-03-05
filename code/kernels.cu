/*!
  \file kernels.cu
  \brief Kernels to be used all over the place
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include <kernels.h>
#include <cstdlib>
#include <stdio.h>
  
//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
#ifdef USE_GPU
__device__ double atomicAdd(double* address, double val)
{
  /*!
    \brief Atomic addition definition for a GPU.
    \param[out] address address of value to add to
    \param[in] val value to be added to address value
  */
  unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
		    __double_as_longlong(val +
					 __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif 

//==========================================================================
arch_global void mapToFace(int M_s, int M_T, int N_s, int* map, scalar* U, scalar* UF){
  /*!
    \brief Map solution to the faces
    \param[in] M_s number of nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] map map from element index to face index
    \param[in] U element solution
    \param[out] UF interface solution
  */
#ifdef USE_CPU
  for(int t = 0; t < M_T; t++){
    for(int j = 0; j < M_s; j++){
      for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int t = blockIdx.x*blkT+threadIdx.z;
  if ( t < M_T){
    int j = threadIdx.x;
    int fc= threadIdx.y;
#endif

	int face;

	for(int d = 0; d < 2; d++){
	  face = ((t*N_F+fc)*2+d)*M_s+j;
	  // Get the interface we will operate on and relate it to the index of U
	  UF[face] = U[map[face]];
	}
	
#ifdef USE_CPU
      }
    }
#endif
  }
}

//==========================================================================
arch_global void mapToElement(int N_s, int N_E, int M_s, int N_N, int* invmap, scalar* Q, scalar* Qtcj){
  /*!
    \brief Map solution from faces to elements
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] M_s number of nodes per interface
    \param[in] N_N number of neighbors per element    
    \param[in] invmap map from face index to element index
    \param[in] Q interface flux
    \param[out] Qtcj interface flux for the elements
  */

#ifdef USE_CPU
  int blk=0;
  for(int e = 0; e < N_E; e++){
    scalar* sol = new scalar[N_F*N_s]; // holds the solution until we update Q
    for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int blk = threadIdx.z; // sol needs room for all the elements in the block
  int e = blockIdx.x*blkE+blk;
  if (e < N_E){
    int i = threadIdx.x;
    int fc= threadIdx.y;
    extern __shared__ scalar sol[];
#endif

    // initialize sol to zero
#ifdef USE_CPU
    for(int i = 0; i < N_s; i++)
#elif USE_GPU
    if(i<N_s)
#endif
      sol[(blk*N_F+fc)*N_s+i] = 0;

    // make sure all threads initialize sol to zero before moving on
#ifdef USE_GPU
    __syncthreads(); 
#endif
    
    // Get the values of Q from the interfaces
#ifdef USE_CPU
    for(int i = 0; i < M_s*N_N; i++){
#elif USE_GPU
    if(i<M_s*N_N){
#endif
      int solidx = invmap[((e*N_F+fc)*M_s*N_N+i)*2+0];
      int qidx   = invmap[((e*N_F+fc)*M_s*N_N+i)*2+1];
#ifdef USE_CPU
      sol[(blk*N_F+fc)*N_s+solidx] += Qtcj[qidx];
#elif USE_GPU
      atomicAdd(&sol[(blk*N_F+fc)*N_s+solidx], Qtcj[qidx]);
#endif
    }

    // make sure all threads complete the atomic add before moving on
#ifdef USE_GPU
    __syncthreads(); 
#endif

    // Push sol to Q
#ifdef USE_CPU
    for(int i = 0; i < N_s; i++)
#elif USE_GPU
    if(i<N_s)
#endif
      Q[(e*N_F+fc)*N_s+i] = sol[(blk*N_F+fc)*N_s+i];

#ifdef USE_CPU
    }// end loop on fields
    delete[] sol;
#endif
  } // end loop on elements
}
   

//==========================================================================
arch_global void redistribute_sf(int N_G, int N_E, scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac){
  /*!
    \brief Take into account the geometry by multiplying with Jacobians
    \param[in] N_G number of gaussian nodes per element
    \param[in] N_E number of elements
    \param[out] sJ s multiplied by Jacobian
    \param[out] fJ f multiplied by Jacobian
    \param[in] s source array
    \param[in] f flux array
    \param[in] J Jacobian
    \param[in] invJac inverse Jacobian (for the fluxes)
  */

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
      for(int fc = 0; fc < N_F; fc++){

#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int g = threadIdx.x;
    int fc= threadIdx.y;
#endif

	scalar j = J[e];
	scalar sol = 0.0; 

	sJ[(e*N_F+fc)*N_G+g] = s[(e*N_F+fc)*N_G+g] * j;
	for(int alpha = 0; alpha < D; alpha++){
	  for(int a = 0; a < D; a++){
	    sol += invJac[((e*N_G+g)*D+alpha)*D+a]*f[((e*N_F+fc)*N_G+g)*D+a] * j;
	  }
	  fJ[((e*N_F+fc)*N_G+g)*D+alpha] = sol;
	  sol = 0;
	}

#ifdef USE_CPU
      }
    }
#endif
  }
}

//==========================================================================
arch_global void redistribute_q(int M_G, int M_T, scalar* qJ, scalar* q, scalar* JF){
  /*!
    \brief Take into account the geometry by multiplying with Jacobians
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_T number of interfaces
    \param[out] qJ q multiplied by Jacobian
    \param[in] q interface flux array
    \param[in] JF Jacobian of the faces
  */

#ifdef USE_CPU
  for(int t = 0; t < M_T; t++){
    for(int g = 0; g < M_G; g++){
      for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int t = blockIdx.x*blkT+threadIdx.z;
  if ( t < M_T){
    int g = threadIdx.x;
    int fc= threadIdx.y;
#endif

      qJ[((t*N_F+fc)*2+0)*M_G+g] = q[((t*N_F+fc)*2+0)*M_G+g]*JF[t*2+0];
      qJ[((t*N_F+fc)*2+1)*M_G+g] = q[((t*N_F+fc)*2+1)*M_G+g]*JF[t*2+1];
  
#ifdef USE_CPU
      }
    }
#endif
  }
}

//==========================================================================
arch_global void addSFQ(int N_s, int N_E, scalar* A, scalar* S, scalar* F, scalar* Q){
  /*!
    \brief Add source, element flux, and interface flux together
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[out] A A=S+F+Q
    \param[in] S source array
    \param[in] F element flux array
    \param[in] Q interface flux
  */

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      for(int fc = 0; fc < N_F; fc++){

#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
    int fc= threadIdx.y;
#endif

	A[(e*N_F+fc)*N_s+i] = S[(e*N_F+fc)*N_s+i] + F[(e*N_F+fc)*N_s+i] + Q[(e*N_F+fc)*N_s+i];

#ifdef USE_CPU
      }
    }
#endif
  }
}


//==========================================================================
//
//  Host C functions
//
//==========================================================================

extern "C" 
void LmapToFace(int M_s, int M_T, int N_s, int* map, scalar* U, scalar* UF){
  /*!
    \brief Host C function to launch mapToFace kernel.
    \param[in] M_s number of nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] map map from element index to face index
    \param[in] U element solution
    \param[out] UF interface solution
    \section Description
    In GPU mode, launches M_T/blkT blocks of M_s x N_F x blkT
    threads. blkT controls the number of interfaces to set on each
    block
  */
#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_s,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  mapToFace arch_args (M_s, M_T, N_s, map, U, UF);
}

extern "C" 
void LmapToElement(int N_s, int N_E, int M_s, int N_N, int* invmap, scalar* Q, scalar* q){
  /*!
    \brief Host C function to launch mapToElement kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] M_s number of nodes per interface
    \param[in] N_N number of neighbors per element
    \param[in] invmap map from face index to element index
    \param[in] Q interface flux
    \param[out] q interface flux for the elements
    \section Description
    In GPU mode, launches N_E/blkE blocks of max(N_s,M_s*N_N) x N_F x
    blkE threads. blkE controls the number of elements to set on each
    block
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(MAX(N_s,M_s*N_N),N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  mapToElement arch_args_array(blkE*N_F*N_s*sizeof(scalar)) (N_s, N_E, M_s, N_N, invmap, Q, q);
}


extern "C" 
void Lredistribute_sf(int N_G, int N_E, scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac){
  /*!
    \brief Host C function to launch redistribute_sf kernel.
    \param[in] N_G number of gaussian nodes per element
    \param[in] N_E number of elements
    \param[out] sJ s multiplied by Jacobian
    \param[out] fJ f multiplied by Jacobian
    \param[in] s source array
    \param[in] f flux array
    \param[in] J Jacobian
    \param[in] invJac inverse Jacobian (for the fluxes)
    \section Description
    In GPU mode, launches N_E/blkE blocks of N_G x N_F x blkE
    threads. blkE controls the number of elements to set on each block
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_G,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  redistribute_sf arch_args (N_G, N_E, sJ, fJ, s, f, J, invJac);
}

extern "C"
void Lredistribute_q(int M_G, int M_T, scalar* qJ, scalar* q, scalar* JF){
  /*!
    \brief Host C function to launch redistribute_q kernel.
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_T number of interfaces
    \param[out] qJ q multiplied by Jacobian
    \param[in] q interface flux array
    \param[in] JF Jacobian of the faces
    \section Description
    In GPU mode, launches M_T/blkT blocks of M_G x N_F x blkT
    threads. blkT controls the number of interfaces to set on each block
  */
#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_G,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  redistribute_q arch_args (M_G, M_T, qJ, q, JF);
}

extern "C" 
void LaddSFQ(int N_s, int N_E, scalar* A, scalar* S, scalar* F, scalar* Q){
  /*!
    \brief Host C function to launch addSFQ kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[out] A A=S+F+Q
    \param[in] S source array
    \param[in] F element flux array
    \param[in] Q interface flux
    \section Description
    In GPU mode, launches N_E/blkE blocks of N_s x N_F x blkE
    threads. blkE controls the number of elements to set on each block
  */

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  addSFQ arch_args (N_s, N_E, A, S, F, Q);
}

