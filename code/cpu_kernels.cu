#include <cpu_kernels.h>
#include <cstdlib>
#include <stdio.h>

//==========================================================================
//
// Limiter prototype function definitions
//
//==========================================================================

arch_device inline int cpu_signum(scalar val){return val>0? 1 : (val<0? -1 : 0);}
arch_device inline scalar cpu_minabs(scalar* c, int n);
arch_device scalar cpu_minmod  (scalar* c, int n); // eq 2.19 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
arch_device scalar cpu_minmod  (scalar a, scalar b); 
arch_device scalar cpu_minmod2 (scalar* c, int n); // eq 2.20 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
#ifdef USE_CPU
scalar cpu_cminmod (scalar* c, int n, scalar eps); // eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
scalar cpu_cminmod2(scalar* c, int n, scalar eps); // eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
#endif
arch_device int cpu_factorial(int n);

arch_device void getTaylorDerivative(int order, int N_s, scalar* T, int mx, int my, int* DxIdx, int* DyIdx, scalar* ddT);
arch_device scalar CellAvg(int N_G, int ioff, scalar* weight, scalar refArea, scalar* powers, int N_s, scalar* T);
  
//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void cpu_equal(int N_s, int N_E, int N_F, scalar* A, scalar* B){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++)
      for(int fc = 0; fc < N_F; fc++)
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
    int fc = threadIdx.y;
#endif

  A[(e*N_F+fc)*N_s+i] = B[(e*N_F+fc)*N_s+i];

  }
}

//==========================================================================
arch_global void cpu_add(int N_s, int N_E, int N_F, scalar* A, scalar* B, scalar c){

  // A = A + c*B
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++)
      for(int fc = 0; fc < N_F; fc++)
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
    int fc = threadIdx.y;
#endif

  A[(e*N_F+fc)*N_s+i] =  A[(e*N_F+fc)*N_s+i] + c*B[(e*N_F+fc)*N_s+i];

  }
}


//==========================================================================
arch_global void cpu_mapToFace_shallow(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){

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

	int idx = -1;
	int face;
	for(int d = 0; d < 2; d++){
	  face= ((t*N_F+fc)*2+d)*M_s+j;
	  idx = map[face];
	  if(idx != -1){
	    UF[face] = U[idx];
	  }
	  else if (idx == -1){
	    if      (fc == 0) UF[((t*N_F+fc)*2+1)*M_s+j] = UF[((t*N_F+fc)*2+0)*M_s+j]; // eta
	    else if (fc == 1) UF[((t*N_F+fc)*2+1)*M_s+j] = -UF[((t*N_F+fc)*2+0)*M_s+j]; // ux
	    else if (fc == 2) UF[((t*N_F+fc)*2+1)*M_s+j] = -UF[((t*N_F+fc)*2+0)*M_s+j]; // uy
	  }
	}
  
#ifdef USE_CPU
      }
    }
#endif
  }
}

//==========================================================================
arch_global void cpu_mapToFace_mhd(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){

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

	int idx = -1;
	int face;

	for(int d = 0; d < 2; d++){
	  face= ((t*N_F+fc)*2+d)*M_s+j;
	  idx = map[face];
	  if(idx != -1){
	    UF[face] = U[idx];
	  }
	}

#ifdef USE_CPU
      }
    }
#endif
  }
}

//==========================================================================
arch_global void cpu_mapToFace(int M_s, int M_T, int N_F, int N_s, int* map, scalar* U, scalar* UF){

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
arch_global void cpu_mapGhostFace(int M_s, int M_ghosts, int N_F, int* ghostInterfaces, scalar* UF){

  /* Communicate the ghost edges in the mesh to and from the other
     partitions

     Just CPU version for now
  */


#ifdef USE_MPI

  MPI_Status status[M_ghosts];
  MPI_Request request[M_ghosts];
  
  MPI_Datatype strided; // make a strided datatype to access UF (only access on side of the interface)
  // from http://stackoverflow.com/questions/15483360/mpi-sending-segments-of-an-array
  // and http://www.mcs.anl.gov/research/projects/mpi/www/www3/MPI_Type_vector.html
  MPI_Type_vector (N_F, M_s, 2*M_s, MPI_SCALAR, &strided);
  MPI_Type_commit(&strided);

  // Allocate a buffer for buffered send
  int bufsize =  M_ghosts * (MPI_BSEND_OVERHEAD + N_F*M_s) + MPI_BSEND_OVERHEAD;
  scalar *buf = new scalar[bufsize];
  MPI_Buffer_attach( buf, bufsize );
  
  int t, dest_source, tag,id;
  //MPI_Comm_rank(MPI_COMM_WORLD,&id);
  //printf("M_ghosts=%i\n",M_ghosts);
  for(int k = 0; k < M_ghosts; k++){
    t           = ghostInterfaces[k*3+0];
    dest_source = ghostInterfaces[k*3+1];
    tag         = ghostInterfaces[k*3+2];

    //printf("CPU id %i is sending/receiving interface %i to proc %i with tag %i\n",id,t,dest_source,tag);

    // VERSION 1: non-blocking send and receive
    // Non-blocking send of my interface
    //MPI_Isend(&UF[t*N_F*2*M_s], 1, strided, dest_source, tag, MPI_COMM_WORLD, &request[2*k+0]);
    // Non-blocking receive of the other side of that interface
    //MPI_Irecv(&UF[t*N_F*2*M_s+M_s], 1, strided, dest_source, tag, MPI_COMM_WORLD, &request[2*k+1]);

    // VERSION 2: buffered send and non-blocking receive
    MPI_Bsend(&UF[t*N_F*2*M_s+M_s], 1, strided, dest_source, tag, MPI_COMM_WORLD);
    MPI_Irecv(&UF[t*N_F*2*M_s+M_s], 1, strided, dest_source, tag, MPI_COMM_WORLD, &request[k]);
  }

  // Wait for all communications to end
  MPI_Waitall(M_ghosts, request, status);
  //printf("Sending/Receiving done by cpu %i\n",id);

  // Detach the buffer
  MPI_Buffer_detach( &buf, &bufsize);
  delete[] buf;
  
  MPI_Type_free(&strided);
#endif
}

//==========================================================================
arch_global void cpu_mapToElement(int N_s, int N_E, int N_F, int M_s, int N_N, int* invmap, scalar* Q, scalar* Qtcj){

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

    // Get the values of Q from the interfaces
#ifdef USE_CPU
    for(int i = 0; i < M_s*N_N; i++){
#elif USE_GPU
    if(i<M_s*N_N){
#endif
      int solidx = invmap[((e*N_F+fc)*M_s*N_N+i)*2+0];
      int qidx   = invmap[((e*N_F+fc)*M_s*N_N+i)*2+1];
      sol[(blk*N_F+fc)*N_s+solidx] += Qtcj[qidx];
    }

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
  } // end loop on elements
#endif
}
   
//   //==========================================================================
// arch_global void cpu_mapToElement(int N_s, int N_E, int N_F, int* invmap, scalar* Q, scalar* Qtcj){

// #ifdef USE_CPU
//   for(int e = 0; e < N_E; e++){
//     for(int i = 0; i < N_s; i++){
//       for(int fc = 0; fc < N_F; fc++){
// #elif USE_GPU
//   int e = blockIdx.x*blkE+threadIdx.z;
//   if (e < N_E){
//     int i = threadIdx.x;
//     int fc= threadIdx.y;
// #endif
// 	int idx = 0;
// 	scalar sol = 0;
// 	for(int k = 0; k < 2; k++){
// 	  idx = invmap[((e*N_F+fc)*2+k)*N_s+i];
// 	  if(idx != -1){
// 	    sol += Qtcj[idx];
// 	  }
// 	}
// 	Q[(e*N_F+fc)*N_s+i] = sol;

// #ifdef USE_CPU
//       }
//     }
// #endif
//   }
// }

//==========================================================================
arch_global void cpu_collocationU(int D, int N_G, int N_s, int N_E, int N_F, scalar* Ug, scalar* dUg, scalar* phi, scalar* dphi, scalar* U){

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

	scalar sol = 0;
	for(int i = 0; i < N_s; i++){
	  sol += phi[i*N_G+g] * U[(e*N_F+fc)*N_s+i];
	}
	Ug[(e*N_F+fc)*N_G+g] = sol;

	sol = 0.0;
	for(int a = 0; a < D; a++){
	  for(int i = 0; i < N_s; i++){
	    sol += dphi[(i*N_G+g)*D+a] * U[(e*N_F+fc)*N_s+i];
	  }
	  dUg[((e*N_F+fc)*N_G+g)*D+a] = sol;
	  sol = 0.0;
	}

#ifdef USE_CPU
      }
    }
#endif
  }
}

//==========================================================================
arch_global void cpu_collocationUF(int M_G, int M_s, int M_T, int N_F, scalar* UgF, scalar* psi, scalar* UF){

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

	scalar sol = 0;
	for(int d = 0; d < 2; d++){
	  for(int j = 0; j < M_s; j++){
	    sol += psi[j*M_G+g] * UF[((t*N_F+fc)*2+d)*M_s+j];
	  }
	  UgF[((t*N_F+fc)*2+d)*M_G+g] = sol;
	  sol = 0.0;
	}

#ifdef USE_CPU
      }
    }
#endif
  }
}



//==========================================================================
arch_global void cpu_redistribute_sf(int D, int N_G, int N_E, int N_F, scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac){

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
arch_global void cpu_gemm_sf(int D, int N_G, int N_s, int N_E, int N_F, scalar* S, scalar* F, scalar* sJ, scalar* fJ, scalar* phi_w, scalar* dphi_w){

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

	scalar sol = 0.0;

	// S = phi_w.transpose() x sJ
	for(int g = 0; g < N_G; g++){
	  sol += phi_w[i*N_G+g] * sJ[(e*N_F+fc)*N_G+g];
	}
	S[(e*N_F+fc)*N_s+i] = sol;
	sol = 0.0;
  
	// F = dphi_w.transpose() x fJ
	sol = 0.0; 
	for(int g = 0; g < N_G; g++){
	  for(int a = 0; a < D; a++){
	    sol += dphi_w[(i*N_G+g)*D+a] * fJ[((e*N_F+fc)*N_G+g)*D+a];
	  }
	}
	F[(e*N_F+fc)*N_s+i] = sol;
	sol = 0.0;

#ifdef USE_CPU
      }
    }
#endif
  }
}

//==========================================================================
arch_global void cpu_redistribute_q(int M_G, int M_T, int N_F, scalar* qJ, scalar* q, scalar* JF){

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
arch_global void cpu_gemm_q(int M_G, int M_s, int M_T, int N_F, scalar* Qtcj, scalar* qJ, scalar* psi_w){

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

	scalar sol = 0.0;
	
	// Qtcj = psi_w.transpose() x qJ
	for(int d = 0; d < 2; d++){
	  for(int g = 0; g < M_G; g++){
	    sol += psi_w[j*M_G+g] * qJ[((t*N_F+fc)*2+d)*M_G+g];
	  }
	  Qtcj[((t*N_F+fc)*2+d)*M_s+j] = sol;
	  sol = 0.0;
	}

#ifdef USE_CPU
      }
    }
#endif
  }
}

//==========================================================================
arch_global void cpu_addSFQ(int N_s, int N_E, int N_F, scalar* A, scalar* S, scalar* F, scalar* Q){

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
arch_global void cpu_zeroVector(int N_s, int N_E, int N_F, scalar* Q){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++)
      for(int fc = 0; fc < N_F; fc++)
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
    int fc= threadIdx.y;
#endif
  
  Q[(e*N_F+fc)*N_s+i] = 0.0;

  }
}

//==========================================================================
arch_global void cpu_hsl(int N_s, int N_E, int N_F, int boundaryMap, scalar* U, scalar* UNew){

  // scalar* c = new scalar[3];
  
  // // Treat boundary conditions
  // if(boundaryMap!=0){ // periodic
  //   for(int fc = 0; fc < N_F; fc++){
  //     for(int m = N_s-1; m > 0; m--){
  // 	// Left side
  // 	c[0] = sqrt((2*m+1)*(2*m+3))*U[(0*N_F+fc)*N_s+m];
  // 	c[1] = U[((0+1)*N_F+fc)*N_s+m-1] - U[( 0     *N_F+fc)*N_s+m-1];
  // 	c[2] = U[( 0   *N_F+fc)*N_s+m-1] - U[((N_E-1)*N_F+fc)*N_s+m-1];
  // 	UNew[(0*N_F+fc)*N_s+m] = 1.0/sqrt((2*m+1)*(2*m+3))*minmod(c,3);
  // 	// Right side
  // 	c[0] = sqrt((2*m+1)*(2*m+3))*U[(0*N_F+fc)*N_s+m];
  // 	c[1] = U[((0    )*N_F+fc)*N_s+m-1] - U[((N_E-1)*N_F+fc)*N_s+m-1];
  // 	c[2] = U[((N_E-1)*N_F+fc)*N_s+m-1] - U[((N_E-2)*N_F+fc)*N_s+m-1];
  // 	UNew[((N_E-1)*N_F+fc)*N_s+m] = 1.0/sqrt((2*m+1)*(2*m+3))*minmod(c,3);
  //     }
  //   }
  // }
  // else if(boundaryMap=0){ // farfield
  //   for(int fc = 0; fc < N_F; fc++){
  //     for(int m = N_s-1; m > 0; m--){
  // 	// Left side
  // 	c[0] = sqrt((2*m+1)*(2*m+3))*U[(0*N_F+fc)*N_s+m];
  // 	c[1] = U[((0+1)*N_F+fc)*N_s+m-1] - U[( 0     *N_F+fc)*N_s+m-1];
  // 	c[2] = 0;
  // 	UNew[(0*N_F+fc)*N_s+m] = 1.0/sqrt((2*m+1)*(2*m+3))*minmod(c,3);
  // 	// Right side
  // 	c[0] = sqrt((2*m+1)*(2*m+3))*U[(0*N_F+fc)*N_s+m];
  // 	c[1] = 0;
  // 	c[2] = U[((N_E-1)*N_F+fc)*N_s+m-1] - U[((N_E-2)*N_F+fc)*N_s+m-1];
  // 	UNew[((N_E-1)*N_F+fc)*N_s+m] = 1.0/sqrt((2*m+1)*(2*m+3))*minmod(c,3);
  //     }
  //   }
  // }
  
  // // Do the other elements
  // for(int e = 1; e < N_E-1; e++){
  //   for(int fc = 0; fc < N_F; fc++){
  //     for(int m = N_s-1; m > 0; m--){
  // 	c[0] = sqrt((2*m+1)*(2*m+3))*U[(e*N_F+fc)*N_s+m];
  // 	c[1] = U[((e+1)*N_F+fc)*N_s+m-1] - U[( e   *N_F+fc)*N_s+m-1];
  // 	c[2] = U[( e   *N_F+fc)*N_s+m-1] - U[((e-1)*N_F+fc)*N_s+m-1];
  // 	UNew[(e*N_F+fc)*N_s+m] = 1.0/sqrt((2*m+1)*(2*m+3))*minmod(c,3);
  // 	if (UNew[(e*N_F+fc)*N_s+m] == U[(e*N_F+fc)*N_s+m]) break;
  //     }
  //   }
  // }

  // delete[] c;
  
}

//==========================================================================
arch_global void cpu_CommunicateGhosts(int N_s, int N_E, int N_F, int N_ghosts, int* ghostElementSend, int* ghostElementRecv, scalar* A){

  /* This function, used by the limiting procedure, communicates the
     elements which are not in my partition. Basically you send the
     elements of A that other partitions need and you receive the
     elements from other partitions that you will need. You store
     these in the back columns of A.*/

#ifdef USE_MPI
  MPI_Status status[2*N_ghosts];
  MPI_Request request[2*N_ghosts];
  int e, dest, source, tag;

  for(int k=0; k<N_ghosts; k++){

    // Send info 
    e   = ghostElementSend[k*3+0]; // local index of A to send
    dest = ghostElementSend[k*3+1]; // destination processor
    tag = ghostElementSend[k*3+2]; // global idx of element
    MPI_Isend(&A[e*N_F*N_s], N_F*N_s, MPI_SCALAR, dest, tag, MPI_COMM_WORLD, &request[2*k+0]);
    
    // Recv info
    e      = ghostElementRecv[k*3+0]; // local index of A to recv (back columns)
    source = ghostElementRecv[k*3+1]; // destination processor
    tag    = ghostElementRecv[k*3+2]; // global idx of element
    MPI_Irecv(&A[e*N_F*N_s], N_F*N_s, MPI_SCALAR, source, tag, MPI_COMM_WORLD, &request[2*k+1]);
  }

  // Wait for communication to end
  MPI_Waitall(2*N_ghosts, request, status);

#endif


}
 
//==========================================================================
arch_global void cpu_hrl1D(int N_s, int N_E, int N_F, int N_G, int N_N, int slicenum, int* neighbors, int offxy, scalar* weight, scalar* V, scalar* A, scalar* Alim){

#ifdef USE_CPU
  int blk = 0;
  for(int e = 0; e < N_E; e++){
    scalar* c = new scalar[2*N_F*slicenum];
    for(int slice = 0; slice < slicenum; slice++){
      for(int fc = 0; fc < N_F; fc++){
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
	      dUL += A[(left *N_F+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*V[j*N_G+g];
	      dUC += A[(e    *N_F+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*V[j*N_G+g];
	      dUR += A[(right*N_F+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*V[j*N_G+g];
	      if(j>=2){
	    	RL += Alim[(e*N_F+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*pow(V[1*N_G+g]-2,j)/(scalar)cpu_factorial(j);
	    	RC += Alim[(e*N_F+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*V[j*N_G+g];
	    	RR += Alim[(e*N_F+fc)*N_s*slicenum+slice*N_s+(j+m-1)]*pow(V[1*N_G+g]+2,j)/(scalar)cpu_factorial(j);
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
	  c[((blk*slicenum+slice)*N_F+fc)*2+0] = 0.5*(avgLC - avgLL);  // 1/dx = 1/2 = 0.5
	  c[((blk*slicenum+slice)*N_F+fc)*2+1] = 0.5*(avgLR - avgLC);

	  Alim[(e*N_F+fc)*N_s*slicenum+slice*N_s+m] = cpu_minmod(&c[((blk*slicenum+slice)*N_F+fc)*2],2); // give it a subset of c (the part you want minmod to act on)
	  //Alim[(e*N_F+fc)*N_s*N_s+slice*N_s+m] = cminmod(c,2,0.01);
	  //or use minmod2(c,2), minmod(c,2,eps), cminmod(c,2,0.01); cminmod2(c,2,eps)
	  if(m==1){Alim[(e*N_F+fc)*N_s*slicenum+slice*N_s+0] = avgLC;}

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
arch_global void cpu_hrl2D(int N_s, int N_E, int N_F, int N_G, int N_N, int D, int order, scalar* XYZCen, scalar* powersXYZG, int* neighbors, int* TaylorDxIdx, int* TaylorDyIdx, scalar* weight, scalar refArea, scalar* A, scalar* Alim){
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
//       	  candA = cpu_minmod(cx,N_N);
//       	  candB = cpu_minmod(cy,N_N);
//       	  if     (k==0)   Alim[(e*N_F+fc)*N_s+m*(m+1)/2+k]   = candA;
//       	  else if(k==m-1){
//       	    Alim[(e*N_F+fc)*N_s+m*(m+1)/2+k]   = cpu_minmod(candA,oldcandA);
//       	    Alim[(e*N_F+fc)*N_s+m*(m+1)/2+k+1] = candB;
//       	  }
//       	  else Alim[(e*N_F+fc)*N_s+m*(m+1)/2+k]= cpu_minmod(candA,oldcandA);
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
arch_global void ChangeBasis(int size1, int size2, int N_E, int N_F, scalar* Transform, scalar* U, scalar* Unew){

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
arch_global void cpu_Prim2Cons(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
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
arch_global void cpu_Cons2Prim(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
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
arch_global void cpu_Half2Cons(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
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
arch_global void cpu_Cons2Half(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
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
void Lcpu_equal(int N_s, int N_E, int N_F, scalar* A, scalar* B){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_equal arch_args (N_s, N_E, N_F, A, B);
}

extern "C" 
void Lcpu_add(int N_s, int N_E, int N_F, scalar* A, scalar* B, scalar c){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif 

  cpu_add arch_args (N_s, N_E, N_F, A, B, c);
}

extern "C" 
void Lcpu_mapToFace_shallow(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){

#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_s,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_mapToFace_shallow arch_args (M_s, M_T, N_F, map, U, UF);
}

extern "C" 
void Lcpu_mapToFace_mhd(int M_s, int M_T, int N_F, int* map, scalar* U, scalar* UF){

#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_s,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_mapToFace_mhd arch_args (M_s, M_T, N_F, map, U, UF);
}

extern "C" 
void Lcpu_mapToFace(int M_s, int M_T, int N_F, int N_s, int* map, scalar* U, scalar* UF){

#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_s,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_mapToFace arch_args (M_s, M_T, N_F, N_s, map, U, UF);
}

extern "C"
void Lcpu_mapGhostFace(int M_s, int M_ghosts, int N_F, int* ghostInterfaces, scalar* UF){

#ifdef USE_GPU
  int div = M_ghosts/blkT;
  int mod = 0;
  if (M_ghosts%blkT != 0) mod = 1;
  dim3 dimBlock(M_s,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_mapGhostFace arch_args (M_s, M_ghosts, N_F, ghostInterfaces, UF);

}

extern "C" 
void Lcpu_mapToElement(int N_s, int N_E, int N_F, int M_s, int N_N, int* invmap, scalar* Q, scalar* q){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(MAX(N_s,M_s*N_N),N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_mapToElement arch_args_array(blkE*N_F*N_s*sizeof(scalar)) (N_s, N_E, N_F, M_s, N_N, invmap, Q, q);
}

extern "C" 
void Lcpu_collocationU(int D, int N_G, int N_s, int N_E, int N_F, scalar* Ug, scalar* dUg, scalar* phi, scalar* dphi, scalar* U){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_G,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_collocationU arch_args (D, N_G, N_s, N_E, N_F, Ug, dUg, phi, dphi, U);
}

extern "C" 
void Lcpu_collocationUF(int M_G, int M_s, int M_T, int N_F, scalar* UgF, scalar* psi, scalar* UF){

#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_G,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_collocationUF arch_args (M_G, M_s, M_T, N_F, UgF, psi, UF);
}

extern "C" 
void Lcpu_redistribute_sf(int D, int N_G, int N_E, int N_F, scalar* sJ, scalar* fJ, scalar* s, scalar* f, scalar* J, scalar* invJac){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_G,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_redistribute_sf arch_args (D, N_G, N_E, N_F, sJ, fJ, s, f, J, invJac);
}

extern "C" 
void Lcpu_gemm_sf(int D, int N_G, int N_s, int N_E, int N_F, scalar* S, scalar* F, scalar* sJ, scalar* fJ, scalar* phi_w, scalar* dphi_w){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_gemm_sf arch_args (D, N_G, N_s, N_E, N_F, S, F, sJ, fJ, phi_w, dphi_w);
}

extern "C"
void Lcpu_redistribute_q(int M_G, int M_T, int N_F, scalar* qJ, scalar* q, scalar* JF){

#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_G,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_redistribute_q arch_args (M_G, M_T, N_F, qJ, q, JF);
}

extern "C" 
void Lcpu_gemm_q(int M_G, int M_s, int M_T, int N_F, scalar* Qtcj, scalar* qJ, scalar* psi_w){

#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_s,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_gemm_q arch_args (M_G, M_s, M_T, N_F, Qtcj, qJ, psi_w);
}

extern "C" 
void Lcpu_addSFQ(int N_s, int N_E, int N_F, scalar* A, scalar* S, scalar* F, scalar* Q){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_addSFQ arch_args (N_s, N_E, N_F, A, S, F, Q);
}

extern "C"
void Lcpu_hsl(int N_s, int N_E, int N_F, int boundaryMap, scalar* U, scalar* UNew){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_hsl arch_args (N_s, N_E, N_F, boundaryMap, U, UNew);
}

extern "C"
void Lcpu_Prim2Cons(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_Prim2Cons arch_args (N_s, N_E, N_F, U, multifluid, passive, model, gamma0);
}

extern "C"
void Lcpu_Cons2Prim(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_Cons2Prim arch_args (N_s, N_E, N_F, U, multifluid, passive, model, gamma0);
}

extern "C"
void Lcpu_Half2Cons(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_Half2Cons arch_args (N_s, N_E, N_F, U, multifluid, passive, model, gamma0);
}

extern "C"
void Lcpu_Cons2Half(int N_s, int N_E, int N_F, scalar* U, bool multifluid, bool passive, int model, scalar gamma0){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_Cons2Half arch_args (N_s, N_E, N_F, U, multifluid, passive, model, gamma0);
}

extern "C"
void Lcpu_CommunicateGhosts(int N_s, int N_E, int N_F, int N_ghosts, int* ghostElementSend, int* ghostElementRecv, scalar* A){

#ifdef USE_GPU
  int div = N_ghosts/blkE;
  int mod = 0;
  if (N_ghosts%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_CommunicateGhosts arch_args (N_s, N_E, N_F, N_ghosts, ghostElementSend, ghostElementRecv, A);

}
  
extern "C"
void Lcpu_hrl1D(int N_s, int N_E, int N_F, int N_G, int N_N, int slicenum, int* neighbors, int offxy, scalar* weight, scalar* V, scalar* A, scalar* Alim){

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(slicenum,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  cpu_hrl1D arch_args_array(blkE*slicenum*N_F*2*sizeof(scalar)) (N_s, N_E, N_F, N_G, N_N, slicenum, neighbors, offxy, weight, V, A, Alim);
}

extern "C"
void Lcpu_hrl2D(int N_s, int N_E, int N_F, int N_G, int N_N, int D, int order, scalar* XYZCen, scalar* powersXYZG, int* neighbors, int* TaylorDxIdx, int* TaylorDyIdx, scalar* weight, scalar refArea, scalar* A, scalar* Alim){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(1,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif
  
  cpu_hrl2D arch_args (N_s, N_E, N_F, N_G, N_N, D, order, XYZCen, powersXYZG, neighbors, TaylorDxIdx, TaylorDyIdx, weight, refArea, A, Alim);
}

extern "C"
void LChangeBasis(int size1, int size2, int N_E, int N_F, scalar* Transform, scalar* U, scalar* Unew){
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(size1,N_F,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  ChangeBasis arch_args (size1, size2, N_E, N_F, Transform, U, Unew);
}

//==========================================================================
//
//  Limiter functions
//
//==========================================================================
arch_device inline scalar cpu_minabs(scalar* c, int n){
  scalar minabs = fabs(c[0]);
  for(int i=1;i<n;i++) if (minabs>fabs(c[i])) minabs = fabs(c[i]);
  return minabs;
}

arch_device scalar cpu_minmod(scalar* c, int n){
  // Generalized minmod function
  // eq 2.19 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  int sign = cpu_signum(c[0]);
  for(int i=1; i<n; i++){
    if (sign!=cpu_signum(c[i])) return 0;
  }
  return sign*cpu_minabs(c,n);
}

arch_device scalar cpu_minmod(scalar a, scalar b){
  // minmod function for 2 arguments
  // eq 2.19 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  int signa = cpu_signum(a);
  if (signa != cpu_signum(b)) return 0;

  scalar fabsa = fabs(a);
  scalar fabsb = fabs(b);
  if (fabsa<fabsb) return signa*fabsa;
  else return signa*fabsb;
}

arch_device scalar cpu_minmod2(scalar* c, int n){
  // eq 2.20 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  scalar min = c[0];
  for(int i=1; i<n; i++) if(fabs(c[i])<fabs(min)) min = c[i];
  return min;
}

#ifdef USE_CPU
scalar cpu_cminmod(scalar* c, int n, scalar eps){
  // eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  // using minmod
  scalar* cc = new scalar[2];
  scalar sum = 0;
  for(int i=0;i<n;i++) sum += c[i];
  cc[0] =(1+eps)*cpu_minmod(c,n);
  cc[1] =(scalar)sum/n;
  scalar m = cpu_minmod(cc,2);
  //delete[] cc;
  return m;    
}

scalar cpu_cminmod2(scalar* c, int n, scalar eps){
  // eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  // using minmod2
  scalar* cc = new scalar[2];
  scalar sum = 0;
  for(int i=0;i<n;i++) sum += c[i];
  cc[0] =(1+eps)*cpu_minmod(c,n);
  cc[1] =(scalar)sum/n;
  scalar m = cpu_minmod2(cc,2);
  delete[] cc;
  return m;    
}
#endif

arch_device int cpu_factorial(int n)
{
#ifdef USE_CPU
  return (n == 1 || n == 0) ? 1 : cpu_factorial(n - 1) * n;
#elif USE_GPU  // no recursion for device less than 2.0
  int f = n;
  while (n>0){
    f*=n;
    n--;
  }
  return f;
#endif 
}

// arch_device void getTaylorDerivative(int order, int N_s, scalar* T, int mx, int my, int* DxIdx, int* DyIdx, scalar* ddT){
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
// }

arch_device scalar CellAvg(int N_G, int ioff, scalar* weight, scalar refArea, scalar* powers, int N_s, scalar* T){
  // Get cell avg of a polynomial of order=order in a cell
  // ioff = 0 for full polynomial
  // ioff = 3 for remainder polynomial
  
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
