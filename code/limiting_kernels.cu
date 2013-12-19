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
  /* Strided copy of array A (length>= numblocks*strideA) to array B (length>= numblocks*strideB)
     numblocks = number of blocks to copy from A to B
     blocklen = number of elements in each block
     strideA = number of elements between start of each block in A
     strideB = number of elements between start of each block in B
     offsetA = number of elements to skip at start of A
     offsetB = number of elements to skip at start of B
     Modeled on MPI_Type_Vector

     You can test with this segment of code:
     scalar* a = new scalar[18];
     scalar* b = new scalar[6];
     for(int i=0; i<18; i++){a[i] = i;printf("%i %f\n",i,a[i]);}
     Lstridedcopy(2,3,9,3,a,b);
     for(int i=0; i<6; i++){printf("%i: %f\n",i,b[i]);}
     delete[] a;
     delete[] b;
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
arch_global void reconstruct_energy(int N_s, int N_E, scalar* rhoeLim, scalar* KLim, scalar* EMono, scalar* ELim){

  /* Reconstruct the energy monomial coefficients using the internal
     and kinetic energy monomial coefficients. Apply a correction to
     the zeroth coefficients so that the method is conservative.*/
     
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
#endif

    int idx, idx0=e*N_s+0;
    
    // Start at idx 1 because we will do the zeroth coefficient separately
    for(int i = 1; i < N_s; i++){
      idx = e*N_s+i;
      ELim[idx] = rhoeLim[idx]+KLim[idx];
    }

    // Correct the zeroth coefficient to conserve energy
    scalar E0 = EMono[idx0];
    for(int k = 2; k<N_s; k+=2){
      idx = e*N_s+k;
      E0 -= 1.0/((scalar)lim_factorial(k+1)) * (ELim[idx]-EMono[idx]);
    }
    ELim[idx0] = E0;
  }
}


//==========================================================================
arch_global void internal_energy(int N_s, int N_E, scalar* p, scalar* g, scalar* rhoe){
  /* Reconstruct the monomial internal energy coefficients using the
     pressure and 1/gamma-1 coefficients so that the pressure remains
     non-oscillatory */
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
#endif

    //printf("==== m = %i\n",i);
    scalar I = 0;
    for(int k=0; k<i+1; k++){
      // could prob do this faster if I brought p and g as a shared array
      //printf("(m,k)=(%i,%i)=%i, m-k=%i, k=%i\n",i,k,binomial_coefficient(i,k),i-k,k);
      I += (scalar)binomial_coefficient(i,k) * p[e*N_s+i-k] * g[e*N_s+k];
    }
    rhoe[e*N_s+i] = I;

#ifdef USE_CPU
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
void Lreconstruct_energy(int N_s, int N_E, scalar* rhoeLim, scalar* KLim, scalar* EMono, scalar* ELim){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(1,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  reconstruct_energy arch_args (N_s, N_E, rhoeLim, KLim, EMono, ELim);
}

extern "C"
void Linternal_energy(int N_s, int N_E, scalar* p, scalar* g, scalar* rhoe){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  internal_energy arch_args (N_s, N_E, p, g, rhoe);

}
 
//==========================================================================
//
//  Internal functions
//
//==========================================================================
arch_device int lim_factorial(int n)
{
#ifdef USE_CPU
  return (n == 1 || n == 0) ? 1 : lim_factorial(n - 1) * n;
#elif USE_GPU  // no recursion for device less than 2.0
  int f = n;
  while (n>0){
    f*=n;
    n--;
  }
  return f;
#endif 
}

arch_device int binomial_coefficient(int n, int k){
  /* Inspired from https://gist.github.com/jeetsukumaran/5392166.
     Does not handle super large numbers (no need really)*/

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

 
