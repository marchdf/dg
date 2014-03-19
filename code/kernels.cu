/*!
  \file kernels.cu
  \brief Kernels to be used all over the place
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
*/
#include <kernels.h>
#include <cstdlib>
#include <stdio.h>

//==========================================================================
//
// Limiter prototype function definitions
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
arch_device int kernel_factorial(int n);

arch_device void getTaylorDerivative(int order, int N_s, scalar* T, int mx, int my, int* DxIdx, int* DyIdx, scalar* ddT);
arch_device scalar CellAvg(int N_G, int ioff, scalar* weight, scalar refArea, scalar* powers, int N_s, scalar* T);
  
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
arch_global void hrl1D(int N_s, int N_E, int N_G, int Nfields, int N_N, int slicenum, int* neighbors, int offxy, scalar* weight, scalar* V, scalar* A, scalar* Alim){
  /*!
    \brief HR limiting function (assumes 1D decomposition)
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] N_G number of gaussian nodes per element
    \param[in] Nfields number of fields to operate on (eg. one field instead of N_F)
    \param[in] N_N number of neighbors per element
    \param[in] slicenum to decompose higher dimensional problem into 1D slices
    \param[in] neighbors array containing an element's neighbors
    \param[in] offxy offset if limiting in x or y
    \param[in] weight gaussian integration weights
    \param[in] V vandermonde matrix evaluated at gaussian nodes
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
	scalar avgdUL = 0, avgdUC=0, avgdUR=0;
	scalar dUL = 0, dUC = 0, dUR = 0;
	scalar RL = 0, RC=0, RR=0;
	scalar avgRL = 0, avgRC=0, avgRR=0;
	scalar avgLL = 0, avgLC=0, avgLR=0;

	int left  = neighbors[e*N_N+offxy+0];
	int right = neighbors[e*N_N+offxy+1];
	//printf("e:%i ; left:%i; right:%i\n",e,neighbors[e*N_N+offxy+0],neighbors[e*N_N+offxy+1]);

	// Check to see if we are at a boundary
	// periodic = 1
	// farfield (i.e. zero-gradient) = 2
	// reflective = 3
	int physical = 0;
	bool L = false; // L = true if I don't have a neighbor on the left
	bool R = false;
	if (left  < 0){physical = -left;  left  = e; L = true;}
	if (right < 0){physical = -right; right = e; R = true;}
	
	// Loop on derivatives
	for(int m = N; m > 0; m--){
	  avgdUL = 0; avgdUC=0; avgdUR=0;

	  // Calculate the derivative average in the cells: left, center,
	  // right calculate the remainder polynomial in our cells and its
	  // two neighbors
	  for(int g = 0; g < N_G; g++){
	    dUL = 0; dUC = 0; dUR = 0;
	    RL  = 0; RC  = 0; RR  = 0;

	    for(int j=0;j<=N-(m-1);j++){
	      dUL += A[(left *Nfields+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*V[j*N_G+g];
	      dUC += A[(e    *Nfields+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*V[j*N_G+g];
	      dUR += A[(right*Nfields+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*V[j*N_G+g];
	      if(j>=2){
	    	RL += Alim[(e*Nfields+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*pow(V[1*N_G+g]-2,j)/(scalar)kernel_factorial(j);
	    	RC += Alim[(e*Nfields+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*V[j*N_G+g];
	    	RR += Alim[(e*Nfields+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*pow(V[1*N_G+g]+2,j)/(scalar)kernel_factorial(j);
	      }// end if
	    }
	    avgdUL += dUL*weight[g];
	    avgdUC += dUC*weight[g];
	    avgdUR += dUR*weight[g];
	    avgRL  += RL*weight[g];
	    avgRC  += RC*weight[g]; 
	    avgRR  += RR*weight[g];
	  }// end integration loop

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
	  
	  avgRL=0; avgRC=0; avgRR=0;
	}// end loop on m

#ifdef USE_CPU
      }
    }
    delete[] c;
#endif
  }
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
arch_global void Prim2Cons(int N_s, int N_E, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
  // Go from primitive to conservative variables

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
      int i = threadIdx.x;
#endif

      scalar rho = U[(e*N_F+0)*N_s+i];
      scalar u   = U[(e*N_F+1)*N_s+i];
      scalar p   = U[(e*N_F+2)*N_s+i];

      if(multifluid){
	scalar gamma= gamma=1.0+1.0/U[(e*N_F+3)*N_s+i];
	U[(e*N_F+1)*N_s+i] = rho*u;
	U[(e*N_F+2)*N_s+i] = p/(gamma-1.0) + 0.5*rho*u*u;
	if      (model==0) U[(e*N_F+3)*N_s+i] = rho/(gamma-1);
      }
      else if(passive){
	scalar phic  = U[(e*N_F+3)*N_s+i];
	U[(e*N_F+1)*N_s+i] = rho*u;
	U[(e*N_F+2)*N_s+i] = p/(gamma0-1) + 0.5*rho*u*u;
	U[(e*N_F+3)*N_s+i] = rho*phic;
      }

#ifdef USE_CPU
    }
#endif
  }
}

//==========================================================================
arch_global void Cons2Prim(int N_s, int N_E, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
  // Go from conservative to primitive variables
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
   int e = blockIdx.x*blkE+threadIdx.z;
   if (e < N_E){
     int i = threadIdx.x;
#endif

      scalar rho  = U[(e*N_F+0)*N_s+i];
      scalar rhou = U[(e*N_F+1)*N_s+i];
      scalar E    = U[(e*N_F+2)*N_s+i];

      if(multifluid){
	scalar gamma=0;
	if      (model==0) gamma=1.0+rho/U[(e*N_F+3)*N_s+i];
	else if (model==1) gamma=1.0+1.0/U[(e*N_F+3)*N_s+i];
	U[(e*N_F+1)*N_s+i] = rhou/rho;
	U[(e*N_F+2)*N_s+i] = (gamma-1.0)*(E - 0.5*rhou*rhou/rho);
	if      (model==0) U[(e*N_F+3)*N_s+i] = 1.0/(gamma-1);
      }
      else if(passive){
	scalar rhophic  = U[(e*N_F+3)*N_s+i];
	U[(e*N_F+1)*N_s+i] = rhou/rho;
	U[(e*N_F+2)*N_s+i] = (gamma0-1)*(E - 0.5*rhou*rhou/rho);
	U[(e*N_F+3)*N_s+i] = rhophic/rho;
      }

#ifdef USE_CPU
    }
#endif
  }
}

//==========================================================================
arch_global void Half2Cons(int N_s, int N_E, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
  // Go from primitive to conservative variables
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
      int i = threadIdx.x;
#endif
      scalar rho = U[(e*N_F+0)*N_s+i];
      scalar u   = U[(e*N_F+1)*N_s+i];
      U[(e*N_F+1)*N_s+i] = rho*u;

#ifdef USE_CPU
    }
#endif
  }
}

//==========================================================================
arch_global void Cons2Half(int N_s, int N_E, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
  // Go from conservative to primitive variables
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
      int i = threadIdx.x;
#endif
      scalar rho  = U[(e*N_F+0)*N_s+i];
      scalar rhou = U[(e*N_F+1)*N_s+i];
      U[(e*N_F+1)*N_s+i] = rhou/rho;
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

extern "C"
void LPrim2Cons(int N_s, int N_E, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  Prim2Cons arch_args (N_s, N_E, U, multifluid, passive, model, gamma0);
}

extern "C"
void LCons2Prim(int N_s, int N_E, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  Cons2Prim arch_args (N_s, N_E, U, multifluid, passive, model, gamma0);
}

extern "C"
void LHalf2Cons(int N_s, int N_E, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  Half2Cons arch_args (N_s, N_E, U, multifluid, passive, model, gamma0);
}

extern "C"
void LCons2Half(int N_s, int N_E, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  Cons2Half arch_args (N_s, N_E, U, multifluid, passive, model, gamma0);
}

extern "C"
void Lhrl1D(int N_s, int N_E, int N_G, int Nfields, int N_N, int slicenum, int* neighbors, int offxy, scalar* weight, scalar* V, scalar* A, scalar* Alim){
  /*!
    \brief Host C function to launch hrl1D kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] N_G number of gaussian nodes per element
    \param[in] Nfields number of fields to operate on (eg. one field instead of N_F)
    \param[in] N_N number of neighbors per element
    \param[in] slicenum to decompose higher dimensional problem into 1D slices
    \param[in] neighbors array containing an element's neighbors
    \param[in] offxy offset if limiting in x or y
    \param[in] weight gaussian integration weights
    \param[in] V vandermonde matrix evaluated at gaussian nodes
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

  hrl1D arch_args_array(blkE*slicenum*Nfields*2*sizeof(scalar)) (N_s, N_E, N_G, Nfields, N_N, slicenum, neighbors, offxy, weight, V, A, Alim);
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
