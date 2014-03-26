/*!
  \file limiting_kernels.cu
  \brief Kernels used by the Limiting class
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup limiting
*/
#include <limiting_kernels.h>
#include <stdio.h>
#include <stdlib.h>

//==========================================================================
//
// Internal prototype function definitions
//
//==========================================================================
arch_device inline int signum(scalar val){/*!Return sign of a scalar*/return val>0? 1 : (val<0? -1 : 0);}
arch_device inline scalar minabs(scalar* c, int n);
arch_device scalar minmod  (scalar* c, int n); // eq 2.19 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
arch_device scalar minmod  (scalar a, scalar b); 
arch_device scalar minmod2 (scalar* c, int n); // eq 2.20 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
#ifdef USE_CPU
scalar cminmod (scalar* c, int n, scalar eps); // eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
scalar cminmod2(scalar* c, int n, scalar eps); // eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
#endif
arch_device inline scalar integrate_monomial_derivative(int k, int n);
arch_device inline scalar integrate_monomial_derivative_bounds(int k, int n,scalar a, scalar b);

arch_device void getTaylorDerivative(int order, int N_s, scalar* T, int mx, int my, int* DxIdx, int* DyIdx, scalar* ddT);
arch_device scalar CellAvg(int N_G, int ioff, scalar* weight, scalar refArea, scalar* powers, int N_s, scalar* T);

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
arch_global void hrl1D(int N_s, int N_E, int Nfields, int N_N, int slicenum, int* neighbors, int offxy, scalar* A, scalar* Alim){
  /*!
    \brief HR limiting function (assumes 1D decomposition)
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] Nfields number of fields to operate on (eg. one field instead of N_F)
    \param[in] N_N number of neighbors per element
    \param[in] slicenum to decompose higher dimensional problem into 1D slices
    \param[in] neighbors array containing an element's neighbors
    \param[in] offxy offset if limiting in x or y
    \param[in] A solution to limit (monomial form)
    \param[out] Alim limited solution (monomial form)
  */
#ifdef USE_CPU
  int blk = 0;
  for(int e = 0; e < N_E; e++){
    scalar* c = new scalar[2*Nfields*slicenum];
    for(int slice = 0; slice < slicenum; slice++){
      for(int fc = 0; fc < Nfields; fc++){
#elif USE_GPU
  int blk = threadIdx.z; // c needs room for all the elements in the block
  int e = blockIdx.x*blkE+blk;
  if (e < N_E){
    int slice= threadIdx.x;
    int fc= threadIdx.y;
    extern __shared__ scalar c[];
#endif  

	int N = N_s-1;
	scalar avgdUL = 0, avgdUC=0, avgdUR=0; scalar integral = 0;
	scalar avgRL = 0, avgRC=0, avgRR=0; scalar alim = 0;
	scalar avgLL = 0, avgLC=0, avgLR=0;
	
	int left  = neighbors[e*N_N+offxy+0];
	int right = neighbors[e*N_N+offxy+1];
	//printf("e:%i ; left:%i; right:%i\n",e,neighbors[e*N_N+offxy+0],neighbors[e*N_N+offxy+1]);

	// Check to see if we are at a boundary
	// periodic = 1
	// farfield (i.e. zero-gradient) = 2
	// reflective = 3
	int physical = 0;
	if (left  < 0){physical = -left;  left  = e;}
	if (right < 0){physical = -right; right = e;}
	
	// Loop on derivatives
	for(int m = N; m > 0; m--){
	  avgdUL = 0; avgdUC=0; avgdUR=0;
	  avgRL = 0; avgRC = 0; avgRR = 0;

	  // Calculate the derivative average in the cells: left,
	  // center, right. Calculate the remainder polynomial in our
	  // cells and its two neighbors
	  for(int n=m-1; n<=N; n++){
	    integral = integrate_monomial_derivative(m-1,n);
	    avgdUL += A[(left *Nfields+fc)*N_s*slicenum+slice*N_s+n]*integral;
	    avgdUC += A[(e    *Nfields+fc)*N_s*slicenum+slice*N_s+n]*integral;
	    avgdUR += A[(right*Nfields+fc)*N_s*slicenum+slice*N_s+n]*integral;
	    if(n>=m+1){
	      alim = Alim[(e*Nfields+fc)*N_s*slicenum+slice*N_s+n];
	      avgRL += alim*integrate_monomial_derivative_bounds(m-1,n,-3,-1);
	      avgRC += alim*integral;
	      avgRR += alim*integrate_monomial_derivative_bounds(m-1,n,1,3);
	    }
	  }
	  
	  // Approximate the average of the linear part
	  avgLL = 0.5*(avgdUL - avgRL); // avg = \frac{1}{2} \int_{-1}^1 U \ud x
	  avgLC = 0.5*(avgdUC - avgRC);
	  avgLR = 0.5*(avgdUR - avgRR);
	
	  // MUSCL approach to get candidate coefficients
	  c[((blk*slicenum+slice)*Nfields+fc)*2+0] = 0.5*(avgLC - avgLL);  // 1/dx = 1/2 = 0.5
	  c[((blk*slicenum+slice)*Nfields+fc)*2+1] = 0.5*(avgLR - avgLC);

	  // farfield (zero-gradient) and reflective BC
	  // Not necessary actually (I think... maybe it is for reflective BC...)
	  if ((physical==2)||(physical==3)){ 
	    Alim[(e*Nfields+fc)*N_s*slicenum+slice*N_s+m] = 0;
	    if(m==1){Alim[(e*Nfields+fc)*N_s*slicenum+slice*N_s+0] = avgLC;}
	  }
	  else if (physical==4){ // gravity field: leave slopes unchanged. This is a problem if we also have shocks...
	    Alim[(e*Nfields+fc)*N_s*slicenum+slice*N_s+m] = A[(e*Nfields+fc)*N_s*slicenum+slice*N_s+m];
	    if(m==1){Alim[(e*Nfields+fc)*N_s*slicenum+slice*N_s+0] = A[(e*Nfields+fc)*N_s*slicenum+slice*N_s+0];}
	  }
	  // Otherwise...
	  else{
	    Alim[(e*Nfields+fc)*N_s*slicenum+slice*N_s+m] = minmod(&c[((blk*slicenum+slice)*Nfields+fc)*2],2); // give it a subset of c (the part you want minmod to act on)
	    //Alim[(e*Nfields+fc)*N_s*N_s+slice*N_s+m] = cminmod(c,2,0.01);
	    //or use minmod2(c,2), minmod(c,2,eps), cminmod(c,2,0.01); cminmod2(c,2,eps)
	    if(m==1){Alim[(e*Nfields+fc)*N_s*slicenum+slice*N_s+0] = avgLC;}
	  }
	}// end loop on m

#ifdef USE_CPU
      }
    }
    delete[] c;
#endif
  }
}

//==========================================================================
arch_global void hri1D(int N_s, int N_E, int N_N, int* neighbors, int N_s1D, int slicenum, int offxy, scalar* Lag2Mono, scalar* Mono2Lag, scalar* U){
  /*!
    \brief HR limiting function for individual elements (assumes 1D decomposition)
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] N_N number of neighbors per element
    \param[in] neighbors array containing an element's neighbors
    \param[in] N_s1D number of nodes in a slice (1D element)
    \param[in] slicenum number of slices (in 2D N_s1D = slicenum)
    \param[in] offxy offset if limiting in x or y
    \param[out] U solution to limit (Lagrange form)
  */

  // Allocations 
  scalar* AL = new scalar[N_s*N_F];
  scalar* AC = new scalar[N_s*N_F];
  scalar* AR = new scalar[N_s*N_F];
  scalar* Alim = new scalar[N_s*N_F];
  scalar sumL=0, sumC=0, sumR=0;

  // Initializations
  scalar avgdUL = 0, avgdUC=0, avgdUR=0; scalar integral = 0;
  scalar avgRL = 0, avgRC=0, avgRR=0; scalar alim = 0;
  scalar avgLL = 0, avgLC=0, avgLR=0;
  scalar c1,c2;
  scalar l2m;
  scalar* uL, *uC, *uR;
  int N = N_s1D-1; // polynomial order
  
  for(int e=0; e<N_E; e++){

    // Neighbors
    int left  = neighbors[e*N_N+offxy+0];
    int right = neighbors[e*N_N+offxy+1];

    // pointer to data of this element
    uL = &U[N_s*left*N_F];
    uC = &U[N_s*e*N_F];
    uR = &U[N_s*right*N_F];
    
    // Check to see if we are at a boundary
    // periodic = 1
    // farfield (i.e. zero-gradient) = 2
    // reflective = 3
    // gravity = 4
    int physical = 0;
    if (left  < 0){physical = -left;  left  = e;}
    if (right < 0){physical = -right; right = e;}

    // Lagrange -> Monomial transformation
    for(int m=0; m<N_s; m++){
      for(int n=0; n<N_F; n++){
	for(int p=0; p<N_s; p++){
	  l2m = Lag2Mono[p*N_s+m];
	  sumL += l2m*uL[n*N_s+p];
	  sumC += l2m*uC[n*N_s+p];
	  sumR += l2m*uR[n*N_s+p];
	}
	AL[n*N_s+m] = sumL; sumL = 0;
	AC[n*N_s+m] = sumC; sumC = 0; Alim[n*N_s+m] = AC[n*N_s+m];
	AR[n*N_s+m] = sumR; sumR = 0;
      }
    }

    if (physical==4){} // gravity field: leave data unchanged. This is a problem if we also have shocks...
    
    // farfield (zero-gradient) and reflective BC
    // Not necessary actually (I think... maybe it is for reflective BC...)
    // set to average in the cell (ie set slopes to 0)
    else if ((physical==2)||(physical==3)){ 
      scalar avgU = 0;
      for(int fc=0; fc<N_F; fc++){                      // loop on fields
	for(int slice = 0; slice < slicenum; slice++){ 	// loop on slices 
	  avgU = 0;
	  // Calculate the cell average and set limited slopes to 0
	  for(int n=0; n<=N; n++){
	    avgU += AC[fc*N_s1D*slicenum+slice*N_s1D+n]*integrate_monomial_derivative(0,n);
	    Alim[fc*N_s1D*slicenum+slice*N_s1D+n] = 0;
	  }

	  // set to cell average
	  Alim[fc*N_s1D*slicenum+slice*N_s1D+0] = 0.5*avgU;
	}
      }
    }

    //Otherwise do the full limiting
    else{
      // loop on fields (in HR: all fields are treated the same way)
      for(int fc=0; fc<N_F; fc++){
	// loop on slices 
	for(int slice = 0; slice < slicenum; slice++){
	
	  // Loop on derivatives
	  for(int m = N; m > 0; m--){
	    avgdUL = 0; avgdUC=0; avgdUR=0;
	    avgRL = 0; avgRC = 0; avgRR = 0;

	    // Calculate the derivative average in the cells: left,
	    // center, right. Calculate the remainder polynomial in our
	    // cells and its two neighbors
	    for(int n=m-1; n<=N; n++){
	      integral = integrate_monomial_derivative(m-1,n);
	      avgdUL += AL[fc*N_s1D*slicenum+slice*N_s1D+n]*integral;
	      avgdUC += AC[fc*N_s1D*slicenum+slice*N_s1D+n]*integral;
	      avgdUR += AR[fc*N_s1D*slicenum+slice*N_s1D+n]*integral;
	      if(n>=m+1){
		alim = Alim[fc*N_s1D*slicenum+slice*N_s1D+n];
		avgRL += alim*integrate_monomial_derivative_bounds(m-1,n,-3,-1);
		avgRC += alim*integral;
		avgRR += alim*integrate_monomial_derivative_bounds(m-1,n,1,3);
	      }
	    }
	  
	    // Approximate the average of the linear part
	    avgLL = 0.5*(avgdUL - avgRL); // avg = \frac{1}{2} \int_{-1}^1 U \ud x
	    avgLC = 0.5*(avgdUC - avgRC);
	    avgLR = 0.5*(avgdUR - avgRR);
	
	    // MUSCL approach to get candidate coefficients
	    c1 = 0.5*(avgLC - avgLL);  // 1/dx = 1/2 = 0.5
	    c2 = 0.5*(avgLR - avgLC);

	    // Limited value
	    Alim[fc*N_s1D*slicenum+slice*N_s1D+m] = minmod(c1,c2); 
	    //if(m==1){Alim[fc*N_s1D*slicenum+slice*N_s1D+0] = avgLC;}
	  }// end loop on m
	  Alim[fc*N_s1D*slicenum+slice*N_s1D+0] = avgLC;
	} // end loop on slices
      } // end loop on fields

    } // end else on physicals
    
    // Monomial -> Lagrange transformation
    for(int m=0; m<N_s; m++){
      for(int n=0; n<N_F; n++){
	for(int p=0; p<N_s; p++){ sumC += Mono2Lag[p*N_s+m]*Alim[n*N_s+p];}
	uC[n*N_s+m] = sumC; sumC = 0;}}
  } // end loop on elements

  
  // Clean up
  delete[] AL;
  delete[] AC;
  delete[] AR;
  delete[] Alim;
  uL = NULL; uC = NULL; uR = NULL;
}

//==========================================================================
arch_global void hrl2D(int N_s, int N_E, int N_G, int N_N, int order, scalar* XYZCen, scalar* powersXYZG, int* neighbors, int* TaylorDxIdx, int* TaylorDyIdx, scalar* weight, scalar refArea, scalar* A, scalar* Alim){
  /*!
    \brief Not used right now. HR limiting function fully 2D (eg. for triangles)
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] N_G number of gaussian nodes per element
    \param[in] N_N number of neighbors per element
    \param[in] order DG order
    \param[in] XYZCen element centroid coordinates
    \param[in] powersXYZG powers of coordinates of gaussian nodes
    \param[in] neighbors array containing an element's neighbors
    \param[in] TaylorDxIdx indices of Taylor polynomial derivatives in x
    \param[in] TaylorDyIdx indices of Taylor polynomial derivatives in y
    \param[in] weight gaussian integration weights
    \param[in] refArea area of reference triangle
    \param[in] A solution to limit (monomial form)
    \param[out] Alim limited solution (monomial form)
  */
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
      int fc= threadIdx.y;
      extern __shared__ scalar c[];
#endif  

//       // Allocate resources locally
//       scalar* localU = new scalar[N_s*(N_N+1)];
//       scalar* neighPowers = new scalar[N_s*N_G*(N_N+1)]; // powers of X in the cell + neighbors
//       scalar* localPowers = new scalar[N_s*N_G*(N_N+1)]; // powers of X to eval cell poly in neighbors
//       scalar* localXYZCen = new scalar[D*(N_N+1)];
//       scalar* dT = new scalar[N_s]; for(int i=0;i<N_s;i++) dT[i]=0;
//       scalar* cx = new scalar[N_N]; // candidates wrt x to limit
//       scalar* cy = new scalar[N_N]; // candidates wrt y to limit
//       scalar* avgdU = new scalar[N_N+1];
//       scalar* avgR  = new scalar[N_N+1];
//       scalar* avgL  = new scalar[N_N+1];
//       scalar candA, candB, oldcandA;
	  
//       // Copy global to local resources (going to have to make these N_F times bigger for GPU)
//       int el = 0; // place holder for the elements
//       for(int nn = 0; nn < N_N+1; nn++){
//       	if(nn==0) el = e; // acting on the target element
//       	else      el = neighbors[e*N_N+nn-1]; // acting on his neighbors
//       	// Copy U
//       	for(int i = 0; i < N_s; i++) localU[nn*N_s+i] = A[(el*N_F+fc)*N_s+i];
//       	// Copy XYZ centroids
//       	//for(int alpha = 0; alpha < D; alpha ++){
//       	scalar cenx = 0, ceny=0;
//       	for(int g = 0; g < N_G; g++){
//       	  cenx += powersXYZG[((e*(N_N+1)+nn)*N_G+g)*N_s+1];
//       	  ceny += powersXYZG[((e*(N_N+1)+nn)*N_G+g)*N_s+2];
//       	}
//       	//localXYZCen[nn*D+alpha] = XYZCen[el*D+alpha];
//       	localXYZCen[nn*D+0] = cenx/N_G;
//       	localXYZCen[nn*D+1] = ceny/N_G;
//       	//}
//       	// Copy powers of XYZ
//       	for(int g = 0; g < N_G; g++){
//       	  for(int i = 0; i < N_s; i++){
//       	    neighPowers[(nn*N_G+g)*N_s+i] = powersXYZG[((el*(N_N+1)+0)*N_G+g)*N_s+i];
//       	    //if((nn==0)&&(e==0)&&(fc==0)) printf("%e ",neighPowers[(nn*N_G+g)*N_s+i]);
//       	    localPowers[(nn*N_G+g)*N_s+i] = powersXYZG[((e*(N_N+1)+nn)*N_G+g)*N_s+i];
//       	  }
//       	  //if((nn==0)&&(e==0)&&(fc==0)) printf("\n");
//       	}
//       }
      
//       // Loop on derivatives
//       for(int m = order; m > 0; m--){
//       	for(int k = 0; k < m; k++){ // loop on combinations of (m-1) order derivatives
//       	  int mx = m-1-k;
//       	  int my = k;
	  
//       	  for(int nn = 0; nn < N_N+1; nn++){ 
//       	    if(nn==0) el = e; // acting on the target element
//       	    else      el = neighbors[e*N_N+nn-1]; // acting on his neighbors

//       	    // Calculate the cell averages of the target polynomial and neighbor polynomials
//       	    for(int i=0;i<N_s;i++) dT[i]=0;
//       	    getTaylorDerivative(order, N_s, &localU[nn*N_s], mx, my, TaylorDxIdx, TaylorDyIdx, dT);
//       	    // if((e==8)&&(fc==0)){
//       	    //   printf("In element %i\n",el);
//       	    //   for(int i=0;i<N_s;i++){printf("    T[%i]=%f ",i,dT[i]);}
//       	    //   printf("\n");
//       	    // }
//       	    avgdU[nn] = CellAvg(N_G, 0, weight, refArea, &neighPowers[nn*N_G*N_s], N_s, dT);
	      
//       	    // Calculate the cell average of the target remainder (with limited coeffs) on element and neighbors
//       	    for(int i=0;i<N_s;i++) dT[i]=0;
//       	    //getTaylorDerivative(order, N_s, &localU[nn*N_s], mx, my, TaylorDxIdx, TaylorDyIdx, dT);
//       	    getTaylorDerivative(order, N_s, &Alim[(e*N_F+fc)*N_s], mx, my, TaylorDxIdx, TaylorDyIdx, dT);
//       	    //if((e==8)&&(fc==0)) for(int i=0;i<N_s;i++) printf("    T[%i]=%f ",i,dT[i]);
//       	    avgR[nn] = CellAvg(N_G, 3, weight, refArea, &localPowers[nn*N_G*N_s], N_s,dT);
      
//       	    // Estimate the cell averages
//       	    avgL[nn] = avgdU[nn] - avgR[nn];

//       	    //if((e==8)&&(fc==0))printf("avgdU=%f avgR=%f avgL=%f\n",avgdU[nn]*(0.5*0.44444),avgR[nn]*(0.5*0.4444444),avgL[nn]*(0.5*0.4444444));
//       	  }

//       	  // store in the coefficient vectors
//       	  for(int nn = 1; nn < N_N+1; nn++){ // loop on element + neighbors
//       	    cx[nn-1] = (avgL[nn] - avgL[0]) / (localXYZCen[nn*D+0]-localXYZCen[0*D+0]);
//       	    cy[nn-1] = (avgL[nn] - avgL[0]) / (localXYZCen[nn*D+1]-localXYZCen[0*D+1]);
//       	    //if((e==8)&&(fc==0)) printf("avgL(%i)=%f, avgL(0)=%f, dL=%f, dx=%f, dy=%f\n", nn,avgL[nn],avgL[0], avgL[nn] - avgL[0], (localXYZCen[nn*D+0]-localXYZCen[0*D+0]), (localXYZCen[nn*D+1]-localXYZCen[0*D+1]));
//       	  }

//       	  //Get the canditate coefficients
//       	  candA = minmod(cx,N_N);
//       	  candB = minmod(cy,N_N);
//       	  if     (k==0)   Alim[(e*N_F+fc)*N_s+m*(m+1)/2+k]   = candA;
//       	  else if(k==m-1){
//       	    Alim[(e*N_F+fc)*N_s+m*(m+1)/2+k]   = minmod(candA,oldcandA);
//       	    Alim[(e*N_F+fc)*N_s+m*(m+1)/2+k+1] = candB;
//       	  }
//       	  else Alim[(e*N_F+fc)*N_s+m*(m+1)/2+k]= minmod(candA,oldcandA);
//       	  oldcandA = candB;
//       	  //if((e==8)&&(fc==0))printf("candA=%f, candB=%f, oldcandA=%f\n",candA,candB,oldcandA);
	  
//       	  // Cell average invariance
//       	  if(m==1) Alim[(e*N_F+fc)*N_s+0] = avgL[0];

//       	} // loop on combinations
//       } // loop on m
      
//       delete[] localU;
//       delete[] neighPowers;
//       delete[] localPowers;
//       delete[] localXYZCen;
//       delete[] dT;
//       delete[] cx;
//       delete[] cy;
//       delete[] avgdU;
//       delete[] avgR;
//       delete[] avgL;
 
#ifdef USE_CPU
    }
#endif
  }
}

//==========================================================================
arch_global void ChangeBasis(int size1, int size2, int N_E, scalar* Transform, scalar* U, scalar* Unew){
  /*!
    \brief Basis transformation (manual). Do not use this. Use BLAS
    \param[in] size1 number of rows of tranform
    \param[in] size2 number of columns of tranform
    \param[in] N_E number of elements
    \param[in] Transform Basis transform matrix (per element)
    \param[in] U solution to transform
    \param[out] Unew transformed solution
  */
    
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < size1; i++){
      for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
    int fc= threadIdx.y;
#endif

  scalar sol = 0.0;
	
  for(int ii = 0; ii < size2; ii++){
    sol += Transform[(e*size1+i)*size2+ii]*U[(e*N_F+fc)*size2+ii];
  }
  Unew[(e*N_F+fc)*size1+i] = sol;
  sol = 0.0;

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

extern "C"
void Lhrl1D(int N_s, int N_E, int Nfields, int N_N, int slicenum, int* neighbors, int offxy, scalar* A, scalar* Alim){
  /*!
    \brief Host C function to launch hrl1D kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] Nfields number of fields to operate on (eg. one field instead of N_F)
    \param[in] N_N number of neighbors per element
    \param[in] slicenum to decompose higher dimensional problem into 1D slices
    \param[in] neighbors array containing an element's neighbors
    \param[in] offxy offset if limiting in x or y
    \param[in] A solution to limit (monomial form)
    \param[out] Alim limited solution (monomial form)
    \section Description
    In GPU mode, launches N_E/blkE blocks of slicenum x Nfields x blkE
    threads. blkE controls the number of elements to set on each block
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(slicenum,Nfields,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  hrl1D arch_args_array(blkE*slicenum*Nfields*2*sizeof(scalar)) (N_s, N_E, Nfields, N_N, slicenum, neighbors, offxy, A, Alim);
}

extern "C" 
void Lhri1D(int N_s, int N_E, int N_N, int* neighbors, int N_s1D, int slicenum, int offxy, scalar* Lag2Mono, scalar* Mono2Lag, scalar* U){
  hri1D(N_s, N_E, N_N, neighbors, N_s1D, slicenum, offxy, Lag2Mono, Mono2Lag, U);
}
 
extern "C"
void Lhrl2D(int N_s, int N_E, int N_G, int N_N, int order, scalar* XYZCen, scalar* powersXYZG, int* neighbors, int* TaylorDxIdx, int* TaylorDyIdx, scalar* weight, scalar refArea, scalar* A, scalar* Alim){
  /*!
    \brief Host C function to launch hrl2D kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] N_G number of gaussian nodes per element
    \param[in] N_N number of neighbors per element
    \param[in] order DG order
    \param[in] XYZCen element centroid coordinates
    \param[in] powersXYZG powers of coordinates of gaussian nodes
    \param[in] neighbors array containing an element's neighbors
    \param[in] TaylorDxIdx indices of Taylor polynomial derivatives in x
    \param[in] TaylorDyIdx indices of Taylor polynomial derivatives in y
    \param[in] weight gaussian integration weights
    \param[in] refArea area of reference triangle
    \param[in] A solution to limit (monomial form)
    \param[out] Alim limited solution (monomial form)
    \section Description
    In GPU mode, launches N_E/blkE blocks of 1 x N_F x blkE
    threads. blkE controls the number of elements to set on each block
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(1,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif
  
  hrl2D arch_args (N_s, N_E, N_G, N_N, order, XYZCen, powersXYZG, neighbors, TaylorDxIdx, TaylorDyIdx, weight, refArea, A, Alim);
}

extern "C"
void LChangeBasis(int size1, int size2, int N_E, scalar* Transform, scalar* U, scalar* Unew){
  /*!
    \brief Host C function to launch ChangeBasis kernel.
    \param[in] size1 number of rows of tranform
    \param[in] size2 number of columns of tranform
    \param[in] N_E number of elements
    \param[in] Transform Basis transform matrix (per element)
    \param[in] U solution to transform
    \param[out] Unew transformed solution
    \section Description
    In GPU mode, launches N_E/blkE blocks of size1 x N_F x blkE
    threads. blkE controls the number of elements to set on each block
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(size1,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  ChangeBasis arch_args (size1, size2, N_E, Transform, U, Unew);
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

 
//==========================================================================
//
//  Limiter functions
//
//==========================================================================
arch_device inline scalar minabs(scalar* c, int n){
  /*!
    \brief Minimum of the absolute value of c
    \param[in] c array to find minabs of
    \param[in] n number of elements in c
    \return minabs of c
  */
  scalar minabs = fabs(c[0]);
  for(int i=1;i<n;i++) if (minabs>fabs(c[i])) minabs = fabs(c[i]);
  return minabs;
}

arch_device scalar minmod(scalar* c, int n){
  /*!
    \brief Generalized minmod function
    \param[in] c array to find minmod of
    \param[in] n number of elements in c
    \return minmod of c
    \section Description
    eq 2.19 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  */
  int sign = signum(c[0]);
  for(int i=1; i<n; i++){
    if (sign!=signum(c[i])) return 0;
  }
  return sign*minabs(c,n);
}

arch_device scalar minmod(scalar a, scalar b){
  /*!
    \brief Minmod function for 2 arguments
    \param[in] a first arg
    \param[in] b second arg
    \return minmod(a,b)
    \section Description
    eq 2.19 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  */
  int signa = signum(a);
  if (signa != signum(b)) return 0;

  scalar fabsa = fabs(a);
  scalar fabsb = fabs(b);
  if (fabsa<fabsb) return signa*fabsa;
  else return signa*fabsb;
}

arch_device scalar minmod2(scalar* c, int n){
  /*!
    \brief Generalized minmod function (alternate)
    \param[in] c array to find minmod of
    \param[in] n number of elements in c
    \return minmod of c
    \section Description
    eq 2.20 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  */
  scalar min = c[0];
  for(int i=1; i<n; i++) if(fabs(c[i])<fabs(min)) min = c[i];
  return min;
}

#ifdef USE_CPU
scalar cminmod(scalar* c, int n, scalar eps){
  /*!
    \brief Generalized minmod function (cpu version)
    \param[in] c array to find minmod of
    \param[in] n number of elements in c
    \param[in] eps smallness comparison
    \return minmod of c
    \section Description
    eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  */
  scalar* cc = new scalar[2];
  scalar sum = 0;
  for(int i=0;i<n;i++) sum += c[i];
  cc[0] =(1+eps)*minmod(c,n);
  cc[1] =(scalar)sum/n;
  scalar m = minmod(cc,2);
  delete[] cc;
  return m;    
}

scalar cminmod2(scalar* c, int n, scalar eps){
  /*!
    \brief Generalized minmod function (alternate, cpu version)
    \param[in] c array to find minmod of
    \param[in] n number of elements in c
    \param[in] eps smallness comparison
    \return minmod of c
    \section Description
    eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  */
  scalar* cc = new scalar[2];
  scalar sum = 0;
  for(int i=0;i<n;i++) sum += c[i];
  cc[0] =(1+eps)*minmod(c,n);
  cc[1] =(scalar)sum/n;
  scalar m = minmod2(cc,2);
  delete[] cc;
  return m;    
}
#endif

arch_device int kernel_factorial(int n)
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
  return 1; // default return for kernel_factorial
}


arch_device inline scalar integrate_monomial_derivative(int k, int n)
{
  /*!
    \brief The integral of the kth derivative of nth order monomial (from -1 to 1)
    \param[in] k kth derivative of the polynomial
    \param[in] n monomial order
    \return \frac{2}{(n-k+1)!} if n-k+1 is odd, 0 otherwise
    Calculates $\int_{-1}^1 \frac{\partial^k}{\partialx^k} \frac{x^n}{n!} \mathrm{d} x$
  */
  int num = n-k+1;
  if (num%2) return 2.0/(scalar)kernel_factorial(num);
  else return 0.0;
}

arch_device inline scalar integrate_monomial_derivative_bounds(int k, int n, scalar a, scalar b)
{
  /*!
    \brief The integral of the kth derivative of nth order monomial.
    \param[in] k kth derivative of the polynomial
    \param[in] n monomial order
    \param[in] a lower integral bound
    \param[in] b upper integral bound
    \return the integral
    Calculates $\int_{a}^{b} \frac{\partial^k}{\partialx^k} \frac{x^n}{n!} \mathrm{d} x$
  */
  int num = n-k+1;
  return (pow(b,num) - pow(a,num))/(scalar)kernel_factorial(num);
}

arch_device void getTaylorDerivative(int order, int N_s, scalar* T, int mx, int my, int* DxIdx, int* DyIdx, scalar* ddT){
  /*!
    \brief Not used right now
  */
//   // Get mx+my order derivative of Taylor polynomial T
//   int offsetx  = 0, offsety = 0;
//   scalar* dT = new scalar[N_s]; for(int i=0; i < N_s; i++) dT[i] = 0;

//   // Find the offsets
//   for(int p = order; p > order-mx; p--) offsetx += (p+1)*(p+2)/2;
//   for(int p = order; p > order-my; p--) offsety += (p+1)*(p+2)/2;
  
//   // Get mx-order order derivative wrt x
//   for(int k = 0; k < (order-mx+1)*(order-mx+2)/2; k++){
//     dT[k] = T[DxIdx[offsetx + k]];
//   }

//   // Get my-order order derivative wrt y
//   for(int k = 0; k < (order-my+1)*(order-my+2)/2; k++){
//     ddT[k] = dT[DyIdx[offsety + k]];
//   }

//   delete[] dT;
}

arch_device scalar CellAvg(int N_G, int ioff, scalar* weight, scalar refArea, scalar* powers, int N_s, scalar* T){
  /*!
    \brief Not used right now. Get cell avg of a polynomial of order=order in a cell
    \return cell average
    \section Description
    ioff = 0 for full polynomial
    ioff = 3 for remainder polynomial
  */  
  scalar I = 0;
  scalar w = 0;
  for(int g = 0; g < N_G; g++){
    w = weight[g];
    for(int i = ioff; i < N_s; i++){
      I += T[i]*w*powers[g*N_s+i];
    }
  }


  // We want the cell average:
  // avg = \frac{1}{|\Omega|} \int_\Omega U(x,y) \ud \Omega
  //     = \frac{1}{|\Omega|} \sum_k w_k U(x_k,y_k) J
  //     = \frac{1}{J |\omega|} \sum_k w_k U(x_k,y_k) J
  //     = \frac{1}{|\omega|} \sum_k w_k U(x_k,y_k)
  return I/refArea; // omega = 1/2 for a triangle
}
