/*!
  \file limiting_kernels.cu
  \brief Kernels used by the Limiting class
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup limiting
*/
#include <limiting_kernels.h>
#include <stdio.h>

//==========================================================================
//
// Internal prototype function definitions
//
//==========================================================================
arch_device int lim_factorial(int n);
arch_device int binomial_coefficient(int n, int k);

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void stridedcopy(int numblocks, int blocklen, int strideA, int strideB, int offsetA, int offsetB, scalar* A, scalar* B){
  /*!
    \brief Strided copy of array A (length>= numblocks*strideA) to array B (length>= numblocks*strideB)
    \param[in] numblocks number of blocks to copy from A to B
    \param[in] blocklen number of elements in each block
    \param[in] strideA number of elements between start of each block in A
    \param[in] strideB number of elements between start of each block in B
    \param[in] offsetA number of elements to skip at start of A
    \param[in] offsetB number of elements to skip at start of B
    \param[in] A source array
    \param[out] B destination array
    \section Description    
    Modeled on MPI_Type_Vector
    
    You can test with this segment of code:
    scalar* a = new scalar[18];
    scalar* b = new scalar[6];
    for(int i=0; i<18; i++){a[i] = i;printf("%i %f\n",i,a[i]);}
    scalar* d_a;
    scalar* d_b;
    cudaMalloc((void**) &d_a,18*sizeof(scalar));
    cudaMalloc((void**) &d_b,6*sizeof(scalar));
    cudaMemcpy(d_a, a, 18*sizeof(scalar), cudaMemcpyHostToDevice);
    Lstridedcopy(2,3,9,3,0,0,d_a,d_b);
    cudaMemcpy(b, d_b, 6*sizeof(scalar), cudaMemcpyDeviceToHost);
    for(int i=0; i<6; i++){printf("%i: %f\n",i,b[i]);}
    delete[] a;
    delete[] b;
    exit(0);
  */

  
#ifdef USE_CPU

  int indA=offsetA,indB=offsetB;

  // Loop on number of blocks
  for(int i = 0; i < numblocks; i++){

    // Copy each block into B
    for(int j = 0; j < blocklen; j++){
      B[indB+j] = A[indA+j];
    }

    indA = indA+strideA;
    indB = indB+strideB;
  }

#elif USE_GPU
  int i = blockIdx.x*blkE+threadIdx.z;
  int indA=offsetA+i*strideA;
  int indB=offsetB+i*strideB;
  if (i < numblocks){
    int j = threadIdx.x;
    B[indB+j] = A[indA+j];
  }    
#endif  
}

//==========================================================================
arch_global void reconstruct_energy(int N_s, int N_E, int slicenum, scalar* rhoeLim, scalar* KLim, scalar* EMono, scalar* ELim){
  /*!
    \brief Reconstruct the energy monomial coefficients
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] slicenum to decompose higher dimensional problem into 1D slices
    \param[in] rhoeLim limited monomial internal energy
    \param[in] KLim limited monomial kinetic energy
    \param[in] EMono monomial total energy
    \param[out] ELim limited monomial total energy
    \section Description
    Reconstruct the energy monomial coefficients using the internal
    and kinetic energy monomial coefficients. Apply a correction tothe
    zeroth coefficients so that the method is conservative.
  */
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int slice = 0; slice < slicenum; slice++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int slice = threadIdx.x;
#endif

    int idx=0;
    int idx0=e*N_s*slicenum+slice*N_s+0;
    
    // Start at idx 1 because we will do the zeroth coefficient separately
    for(int i = 1; i < N_s; i++){
      idx = e*N_s*slicenum+slice*N_s+i;
      ELim[idx] = rhoeLim[idx]+KLim[idx];
    }

    // Correct the zeroth coefficient to conserve energy
    scalar E0 = EMono[idx0];
    for(int k = 2; k<N_s; k+=2){
      idx = e*N_s*slicenum+slice*N_s+k;
      E0 -= 1.0/((scalar)lim_factorial(k+1)) * (ELim[idx]-EMono[idx]);
    }
    ELim[idx0] = E0;

#ifdef USE_CPU
  }// for on slice
#endif
  }
}


//==========================================================================
arch_global void internal_energy_multifluid(int N_s, int N_E, int slicenum, scalar* p, scalar* g, scalar* rhoe){
  /*!
    \brief Reconstruct the energy monomial coefficients
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] slicenum to decompose higher dimensional problem into 1D slices
    \param[in] p monomial pressure solution
    \param[in] g monomial 1/(gamma-1) solution
    \param[out] rhoe monomial internal energy
    \section Description
    Reconstruct the monomial internal energy coefficients using the
    pressure and 1/gamma-1 coefficients so that the pressure remains
    non-oscillatory
  */
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int slice = 0; slice < slicenum; slice++){
      for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
    int slice = threadIdx.y;
#endif

    //printf("==== m = %i\n",i);
    scalar I = 0;
    for(int k=0; k<i+1; k++){
      // could prob do this faster if I brought p and g as a shared array
      //printf("(m,k)=(%i,%i)=%i, m-k=%i, k=%i\n",i,k,binomial_coefficient(i,k),i-k,k);
      I += (scalar)binomial_coefficient(i,k) * p[e*N_s*slicenum+slice*N_s+i-k] * g[e*N_s*slicenum+slice*N_s+k];
    }
    rhoe[e*N_s*slicenum+slice*N_s+i] = I;

#ifdef USE_CPU
    }
  }
#endif
  }
}

//==========================================================================
arch_global void internal_energy_stiffened(int N_s, int N_E, int slicenum, scalar* p, scalar* g, scalar* b, scalar* rhoe){
  /*!
    \brief Reconstruct the energy monomial coefficients
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] slicenum to decompose higher dimensional problem into 1D slices
    \param[in] p monomial pressure solution
    \param[in] g monomial 1/(gamma-1) solution
    \param[in] b monomial gamma*pinf/(gamma-1) solution
    \param[out] rhoe monomial internal energy
    \section Description
    Reconstruct the monomial internal energy coefficients using the
    pressure, 1/gamma-1, and gamma*pinf/(gamma-1) coefficients so that
    the pressure remains non-oscillatory
  */
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int slice = 0; slice < slicenum; slice++){
      for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
    int slice = threadIdx.y;
#endif

    //printf("==== m = %i\n",i);
    scalar I = 0;
    for(int k=0; k<i+1; k++){
      // could prob do this faster if I brought p and g as a shared array
      //printf("(m,k)=(%i,%i)=%i, m-k=%i, k=%i\n",i,k,binomial_coefficient(i,k),i-k,k);
      I += (scalar)binomial_coefficient(i,k) * p[e*N_s*slicenum+slice*N_s+i-k] * g[e*N_s*slicenum+slice*N_s+k];
    }
    rhoe[e*N_s*slicenum+slice*N_s+i] = I + b[e*N_s*slicenum+slice*N_s+i];

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
void Lstridedcopy(int numblocks, int blocklen, int strideA, int strideB, int offsetA, int offsetB, scalar* A, scalar* B){
  /*!
    \brief Host C function to lauch stridedcopy kernel.
    \param[in] numblocks number of blocks to copy from A to B
    \param[in] blocklen number of elements in each block
    \param[in] strideA number of elements between start of each block in A
    \param[in] strideB number of elements between start of each block in B
    \param[in] offsetA number of elements to skip at start of A
    \param[in] offsetB number of elements to skip at start of B
    \param[in] A source array
    \param[out] B destination array
    \section Description
    In GPU mode, launches numblocks/blkE blocks of blocklen x 1 x blkE
    threads. blkE controls the number of elements to set on each block
  */
#ifdef USE_GPU
  int div = numblocks/blkE;
  int mod = 0;
  if (numblocks%blkE != 0) mod = 1;
  dim3 dimBlock(blocklen,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  stridedcopy arch_args (numblocks, blocklen, strideA, strideB, offsetA, offsetB, A, B);
};

extern "C"
void Lreconstruct_energy(int N_s, int N_E, int slicenum, scalar* rhoeLim, scalar* KLim, scalar* EMono, scalar* ELim){
  /*!
    \brief Host C function to lauch reconstruct_energy kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] slicenum to decompose higher dimensional problem into 1D slices
    \param[in] rhoeLim limited monomial internal energy
    \param[in] KLim limited monomial kinetic energy
    \param[in] EMono monomial total energy
    \param[out] ELim limited monomial total energy
    \section Description
    In GPU mode, launches N_E/blkE blocks of slicenum x 1 x blkE
    threads. blkE controls the number of elements to set on each block
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(slicenum,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  reconstruct_energy arch_args (N_s, N_E, slicenum, rhoeLim, KLim, EMono, ELim);
}

extern "C"
void Linternal_energy_multifluid(int N_s, int N_E, int slicenum, scalar* p, scalar* g, scalar* rhoe){
  /*!
    \brief Host C function to lauch internal_energy_multifluid kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] slicenum to decompose higher dimensional problem into 1D slices
    \param[in] p monomial pressure solution
    \param[in] g monomial 1/(gamma-1) solution
    \param[out] rhoe monomial internal energy
    \section Description
    In GPU mode, launches N_E/blkE blocks of N_s x slicenum x blkE
    threads. blkE controls the number of elements to set on each block
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,slicenum,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  internal_energy_multifluid arch_args (N_s, N_E, slicenum, p, g, rhoe);

}

extern "C"
void Linternal_energy_stiffened(int N_s, int N_E, int slicenum, scalar* p, scalar* g, scalar* b, scalar* rhoe){
  /*!
    \brief Host C function to lauch internal_energy_multifluid kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] slicenum to decompose higher dimensional problem into 1D slices
    \param[in] p monomial pressure solution
    \param[in] g monomial 1/(gamma-1) solution
    \param[in] b monomial gamma*pinf/(gamma-1) solution
    \param[out] rhoe monomial internal energy
    \section Description
    In GPU mode, launches N_E/blkE blocks of N_s x slicenum x blkE
    threads. blkE controls the number of elements to set on each block
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,slicenum,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  internal_energy_stiffened arch_args (N_s, N_E, slicenum, p, g, b, rhoe);

}

//==========================================================================
//
//  Internal functions
//
//==========================================================================
arch_device int lim_factorial(int n)
{
  /*!
    \brief Factorial function
    \param[in] n get factorial of this number
    \return factorial of n
  */
  if     (n== 0) return 1;
  else if(n== 1) return 1;
  else if(n== 2) return 2;
  else if(n== 3) return 6;
  else if(n== 4) return 24;
  else if(n== 5) return 120;
  else if(n== 6) return 720;
  else if(n== 7) return 5040;
  else if(n== 8) return 40320;
  else if(n== 9) return 362880;
  else if(n==10) return 3628800;
  else if(n==11) return 39916800;
  else if(n==12) return 479001600;
  return 1; // default return for factorial
}

arch_device int binomial_coefficient(int n, int k){
  /*!
    \brief Binomial coefficient function
    \param[in] n
    \param[in] k
    \return C(n,k)
    \section Description
    Inspired from https://gist.github.com/jeetsukumaran/5392166.
    Does not handle super large numbers (no need really)
  */

  if (0 == k || n == k) {
    return 1;
  }
  if (k > n) {
    return 0;
  }
  if (k > (n - k)) {
    k = n - k;
  }
  if (1 == k) {
    return n;
  }
  int b = 1;
  for (int i = 1; i <= k; ++i) {
    b *= (n - (k - i));
    if (b < 0) return -1; /* Overflow */
    b /= i;
  }
  return b;
}

 
