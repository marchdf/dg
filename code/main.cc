#include <string>
#include <iomanip>
#include <stdio.h>
#include "stdlib.h"
#include <cutil.h>
#include <cutil_inline.h>
#include <cublas.h>
#include <time.h>
#include "fullMatrix.h"
#include "polynomialBasis.h"
#include "quadratures/Gauss.h"
#include "GmshDefines.h"
#include "simpleMesh.h"
#include <blas_stuff.h>
#include <scalar_def.h>
#include <dg_functions.h>
#include <gpu_kernels.h>
#include <cpu_kernels.h>
#include <deck.h>
#include <init_cond.h>
#include <print_sol.h>


//
// Function prototypes
//
void copyMatrixToPointer(fullMatrix<scalar> &A, scalar* h_A);
void copyMatrixToPointer(const fullMatrix<scalar> &A, scalar* h_A);
void blasScopy(int N, float* x, int INCX, float* y, int INCY){
  F77NAME(scopy)(&N, x, &INCX, y, &INCY);
}
void blasSaxpy(int M, float alpha, float* x, int INCX, float* y, int INCY){
  F77NAME(saxpy)(&M, &alpha, x ,&INCX, y, &INCY);
}
void blasSgemm(char* or1, char* or2, int M , int N, int K, float alpha, float* A, int LDA, float* B, int LDB, float beta, float* C, int LDC){
  F77NAME(sgemm)(or1, or2, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
void blasDcopy(int N, double* x, int INCX, double* y, int INCY){
  F77NAME(dcopy)(&N, x, &INCX, y, &INCY);
}
void blasDaxpy(int M, double alpha, double* x, int INCX, double* y, int INCY){
  F77NAME(daxpy)(&M, &alpha, x ,&INCX, y, &INCY);
}
void blasDgemm(char* or1, char* or2, int M , int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double beta, double* C, int LDC){
  F77NAME(dgemm)(or1, or2, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
void get_element_types(const int order, int &msh_lin){
  if      (order==0)  { msh_lin = MSH_LIN_2;  }
  else if (order==1)  { msh_lin = MSH_LIN_2;  }
  else if (order==2)  { msh_lin = MSH_LIN_3;  }
  else if (order==3)  { msh_lin = MSH_LIN_4;  }
  else if (order==4)  { msh_lin = MSH_LIN_5;  }
  else if (order==5)  { msh_lin = MSH_LIN_6;  }
  else if (order==6)  { msh_lin = MSH_LIN_7;  }
  else if (order==7)  { msh_lin = MSH_LIN_8;  }
  else if (order==8)  { msh_lin = MSH_LIN_9;  }
  else if (order==9)  { msh_lin = MSH_LIN_10; }
  else if (order==10) { msh_lin = MSH_LIN_11; }
  else {printf("Invalid order number.");}
}
void makeZero(scalar* A, int size){
  for(int k=0; k < size; k++) A[k] = 0.0;
}
void average_cell_p0(const int N_s, const int N_E, const int N_F, fullMatrix<scalar> &U);

int main (int argc, char **argv)
{


  
  ////////////////////////////////////////////////////////////////////////////
  //
  // Read arguments from deck
  //
  ////////////////////////////////////////////////////////////////////////////   
  if (argc!=3)
    printf("No input deck given. Defaulting to deck.inp\n");
  std::string deckfile = "deck.inp";
  for (int i=2;i<argc;i++){
    std::string argType = argv[i-1];
    if (argType== "-d") deckfile = argv[i];
  }
  deck inputs;  
  inputs.readDeck(deckfile.c_str());
  // Choose the processor type (cpu or gpu) cpu by default
  bool cpu = true;
  if (inputs.getProc() == "cpu")      {cpu = true; printf("Using CPU\n");}
  else if (inputs.getProc() == "gpu") {cpu = false;printf("Using GPU\n");}
  else{ printf("Invalid processor type. Correct the deck. Defaulting to cpu.\n");}
  
  bool debug = inputs.getDebug();
  bool blas = inputs.getBlas(); if (blas==1) printf("Using BLAS\n");
  int order = inputs.getOrder();
  bool order0 = false; if (order==0) {order0 = true; order = 1;}

  // Get the flux
  int  flux;
  if      (inputs.getFlux()=="llf") flux = 0;
  else if (inputs.getFlux()=="ncf") flux = 1;
  else if (inputs.getFlux()=="roe") flux = 2;
  else{ printf("Invalid flux setup. Correct the deck.\n");}
  
  std::string fileName = inputs.getMeshfile();
  
  // Setup the problem type
  bool multifluid = false;
  bool passive = false;
  if (inputs.getProblem() == "multifluid")   multifluid = true;
  if (inputs.getProblem() == "passive")      passive = true;
  else{ printf("Invalid problem setup. Correct the deck.\n");}

  // Setup the model
  int model = 0;
  if      (inputs.getModel() == "gammamod")   model = 0;
  else if (inputs.getModel() == "invgamma")   model = 1;
  else if (passive) model = 0;
  else{ printf("Invalid model setup. Correct the deck.\n");}
    
  // Setup the initial condition type
  bool simplew = false;
  bool sodtube = false;
  bool contact = false;
  bool matfrnt = false;
  bool sinegam = false;
  bool expogam = false;
  bool sinephi = false;
  if      (inputs.getInitialCondition()=="simplew") simplew = true;
  else if (inputs.getInitialCondition()=="sodtube") sodtube = true;
  else if (inputs.getInitialCondition()=="contact") contact = true;
  else if (inputs.getInitialCondition()=="matfrnt") matfrnt = true;
  else if (inputs.getInitialCondition()=="sinegam") sinegam = true;
  else if (inputs.getInitialCondition()=="expogam") expogam = true;
  else if (inputs.getInitialCondition()=="sinephi") sinephi = true;
  else{printf("Invalid initial condition setup. Correct the deck.\n");}

  // setup the boundary condition type
  bool periodic = false;
  bool farfield = false;
  if      (inputs.getBoundaryCondition()=="periodic") periodic = true;
  else if (inputs.getBoundaryCondition()=="farfield") farfield = true;
  else{ printf("Invalid boundary condition setup. Correct the deck.\n");}    


  //==========================================================================
  //
  //   CPU calculations
  //
  //==========================================================================


  ////////////////////////////////////////////////////////////////////////////
  //
  // Load a 2D mesh, node coordinates, elements, interfaces, normals
  //
  ////////////////////////////////////////////////////////////////////////////   
  simpleMesh m;
  m.load(fileName.c_str());
  int msh_lin;
  get_element_types(order, msh_lin);
  
  // Get the nodes, elements, interfaces, normals
  const fullMatrix<double> &nodes = m.getNodes();
  const std::vector<simpleElement> &elements = m.getElements(msh_lin);

  ////////////////////////////////////////////////////////////////////////////   
  //
  // Generer les fonctions de formes, get integration points, weights
  //
  ////////////////////////////////////////////////////////////////////////////   

  const polynomialBasis *basis  = polynomialBases::find (msh_lin);  // for the element
  fullMatrix<double> points, weight;
  gaussIntegration::getLine(order*2+1, points, weight);


  //////////////////////////////////////////////////////////////////////////   
  //
  // Define some numbers for clarity
  //
  //////////////////////////////////////////////////////////////////////////   

  int N_s = elements[0].getNbNodes(); // number of nodes on an element          (i index)
  int M_s = 1;                        // number of nodes on a face              (j index)
  int N_T = basis->numFaces;          // number of faces per element            
  int N_E = elements.size();          // number of elements                     (e index)
  int M_T = N_E+1;                    // number of faces                        (t index)
  int N   = N_E*N_s;                  // number of dof for a DG
  int N_G = points.size1();           // number of integration points           (g index)
  int M_G = 1;                        // number of integration points on a face (g index)
  int D   = 1;                        // number of coordinates in elem ref space(alpha index)
  int DF  = 1;                        // number of coordinates in face ref space(alpha index)
  int N_F = inputs.getNumFields();    // number of unknown fields               (h, u_x, u_y with fc index)

  if(order0) printf("order %i\n",0); else printf("order %i\n",order);
  printf("N_s %i\n",N_s);
  printf("M_s %i\n",M_s);
  printf("N_T %i\n",N_T);
  printf("N_E %i\n",N_E);
  printf("M_T %i\n",M_T);
  printf("N_G %i\n",N_G);
  printf("M_G %i\n",M_G);

 
  //////////////////////////////////////////////////////////////////////////   
  //
  // Calcul de la valeur des fonctions de formes aux points d'integration
  //
  //////////////////////////////////////////////////////////////////////////

  // Elements
  fullMatrix<scalar> phi (N_G,N_s); 
  fullMatrix<double> phiD (N_G,N_s); 
  fullMatrix<scalar> dphi(N_G*D,N_s);
  basis->f (points, phiD);
  for(int g = 0; g < N_G; g++){
    for(int i = 0; i < N_s; i++){
      phi(g,i) = (scalar)phiD(g,i);
    }	  
  }    
  double grads[N_s][3];  
  for(int g = 0; g < N_G; g++){
    basis->df(points(g,0),points(g,1),points(g,2),grads);
    for(int alpha = 0; alpha < D; alpha ++){
      for(int i = 0; i < N_s; i++){
  	dphi(g*D+alpha,i) = (scalar)grads[i][alpha];  //see paper for indexing p.6
      }	  
    }    
  }

  // scalar sum = 0.0;
  // for(int g = 0; g < N_G; g++){
  //   sum = 0.0;
  //   for(int i = 0; i < N_s; i++){
  //     sum+=dphi(g,i);  //see paper for indexing p.6
  //   }
  //   printf("g:%i -> sum dphi=%12.11f\n", g, sum); 
  // }
  

  //////////////////////////////////////////////////////////////////////////   
  //
  // Multiply des fonctions de formes et derivees avec les poids
  //
  //////////////////////////////////////////////////////////////////////////

  // Elements
  fullMatrix<scalar> phi_w (N_G,N_s); 
  fullMatrix<scalar> dphi_w(N_G*D,N_s);
  for(int g = 0; g < N_G; g++){
    for(int i = 0; i < N_s; i++){
      phi_w(g,i) = (scalar)phi(g,i) * weight(g,0);
      for(int alpha = 0; alpha < D; alpha++){
  	dphi_w(g*D+alpha,i) = (scalar)dphi(g*D+alpha,i)*weight(g,0);
      }
    }
  }
  
  
  //////////////////////////////////////////////////////////////////////////   
  //
  // Get XYZ coordinates at integration points of each elements and interface
  //
  //////////////////////////////////////////////////////////////////////////
  
  // Elements
  fullMatrix<scalar> XYZNodes (N_s, N_E*D);
  fullMatrix<scalar> XYZG (N_G, N_E*D);
  for (int e = 0; e < N_E; e++) {
    const simpleElement &el = elements[e];
    for (int i = 0; i < N_s; i++) {
      int inode = el.getNode (i);
      for(int alpha = 0; alpha < D; alpha++){
  	XYZNodes (i, e*D+alpha) = (scalar)nodes (alpha, inode);
      }
    }
  }
  XYZG.gemm (phi, XYZNodes);

  //////////////////////////////////////////////////////////////////////////   
  //
  // Build the boundary map
  //
  //////////////////////////////////////////////////////////////////////////

  int boundaryMap;
  if (periodic) boundaryMap = M_T-1;
  if (farfield) boundaryMap = 0; 

  //////////////////////////////////////////////////////////////////////////   
  //
  // Initialize the unknowns
  //
  //////////////////////////////////////////////////////////////////////////
  
  fullMatrix<scalar> U(N_s, N_E*N_F);
  fullMatrix<scalar> Us(N_s, N_E*N_F);
  fullMatrix<scalar> Ustar(N_s, N_E*N_F);
  scalar gamma0 = 0;
  if(multifluid){
    if     (simplew) init_dg_simplew_multifluid(N_s, N_E, N_F, D, model, XYZNodes, U);
    else if(sodtube) init_dg_sodtube_multifluid(N_s, N_E, N_F, D, model, XYZNodes, U);
    else if(contact) init_dg_contact_multifluid(N_s, N_E, N_F, D, model, XYZNodes, U);
    else if(matfrnt) init_dg_matfrnt_multifluid(N_s, N_E, N_F, D, model, XYZNodes, U);
    else if(sinegam) init_dg_sinegam_multifluid(N_s, N_E, N_F, D, model, XYZNodes, U);
    else if(expogam) init_dg_expogam_multifluid(N_s, N_E, N_F, D, model, XYZNodes, U);
  }
  else if (passive){
    if (sinephi) init_dg_sinephi_passive(N_s, N_E, N_F, D, gamma0, XYZNodes, U);
  }

  if (order0) average_cell_p0(N_s, N_E, N_F, U);
  
  // print the initial condition to the file
  printf("Initial condition written to output file.\n");
  if(multifluid)print_dg_multifluid(N_s, N_E, N_F, model, U, m, msh_lin, 0, 0, 0,-1);
  if(passive)   print_dg_passive(N_s, N_E, N_F, gamma0, U, m, msh_lin, 0, 0, 0,-1);

  
  //////////////////////////////////////////////////////////////////////////   
  //
  // Calculate the jacobian matrix for each integration point
  //
  //////////////////////////////////////////////////////////////////////////
  
  // Elements
  fullMatrix<scalar> Jac(N_E*D,N_G*D);
  fullMatrix<scalar> invJac(N_E*D,N_G*D); // Inverse Jacobian matrix: dxi/dx
  fullMatrix<scalar> J(N_E,1);            // determinant of the Jacobian
  fullMatrix<scalar> invJ(N_E,1);         // determinant of the inverse Jacobian
  dg_jacobians_elements(N_G, N_E, D, XYZNodes, dphi, Jac, invJac, J, invJ);

  //////////////////////////////////////////////////////////////////////////   
  // 
  // Calculate the inverse mass matrices
  //
  //////////////////////////////////////////////////////////////////////////
  
  scalar* h_Minv = new scalar[N_s*N_s*N_E];
  dg_inverse_mass_matrix(order, msh_lin, inputs.getElemType(), N_s, N_E, D, XYZNodes, h_Minv);
  
  //==========================================================================
  //
  //   GPU calculations
  //
  //==========================================================================

  // Choose the device
  CUDA_SAFE_CALL(cudaSetDevice(0));

  // Use cublas or not
  cublasStatus status;
  status = cublasInit();

  //////////////////////////////////////////////////////////////////////////   
  //
  // Send stuff to the GPU
  //
  //////////////////////////////////////////////////////////////////////////   

  //
  //  Initialize some stuff on the host. We need to transform
  //  the data in fullMatrix to a pointer form to transfer to GPU
  //  note: it's got to be column major sorted
  //
  scalar* h_phi     = new scalar[N_G*N_s];          makeZero(h_phi,N_G*N_s);
  scalar* h_phi_w   = new scalar[N_G*N_s];          makeZero(h_phi_w,N_G*N_s);          
  scalar* h_dphi    = new scalar[D*N_G*N_s];	    makeZero(h_dphi,D*N_G*N_s);	 
  scalar* h_dphi_w  = new scalar[D*N_G*N_s];	    makeZero(h_dphi_w,D*N_G*N_s);	 
  scalar* h_J       = new scalar[N_E];              makeZero(h_J,N_E);                                 // not same as J!!
  scalar* h_invJac  = new scalar[N_G*D*N_E*D];      makeZero(h_invJac,N_G*D*N_E*D);                    // not same as invJac!!
  scalar* h_Us      = new scalar[N_s*N_E*N_F];	    makeZero(h_Us,N_s*N_E*N_F);	 
  scalar* h_Ustar   = new scalar[N_s*N_E*N_F];	    makeZero(h_Ustar,N_s*N_E*N_F);
  scalar* h_DU      = new scalar[N_s*N_E*N_F];	    makeZero(h_DU,N_s*N_E*N_F);	 
  scalar* h_U       = new scalar[N_s*N_E*N_F];	    makeZero(h_U,N_s*N_E*N_F);	 
  scalar* h_UF      = new scalar[2*N_F*M_s*M_T];    makeZero(h_UF,2*N_F*M_s*M_T); 
  scalar* h_Uinteg  = new scalar[N_F*N_G*N_E];	    makeZero(h_Uinteg,N_F*N_G*N_E);	 
  scalar* h_dUinteg = new scalar[D*N_G*N_E*N_F];    makeZero(h_dUinteg,D*N_G*N_E*N_F); 
  scalar* h_UintegF = new scalar[2*N_F*M_G*M_T];    makeZero(h_UintegF,2*N_F*M_G*M_T); 
  scalar* h_s       = new scalar[N_G*N_E*N_F];	    makeZero(h_s,N_G*N_E*N_F);	 
  scalar* h_sJ      = new scalar[N_G*N_E*N_F];	    makeZero(h_sJ,N_G*N_E*N_F);	 
  scalar* h_S       = new scalar[N_s*N_E*N_F];	    makeZero(h_S,N_s*N_E*N_F);	 
  scalar* h_f       = new scalar[D*N_F*N_G*N_E];    makeZero(h_f,D*N_F*N_G*N_E); 
  scalar* h_fJ      = new scalar[D*N_G*N_E*N_F];    makeZero(h_fJ,D*N_G*N_E*N_F); 
  scalar* h_F       = new scalar[N_s*N_E*N_F];	    makeZero(h_F,N_s*N_E*N_F);	 
  scalar* h_q       = new scalar[M_G*M_T*N_F*2];    makeZero(h_q,M_G*M_T*N_F*2); 
  scalar* h_Q       = new scalar[N_s*N_E*N_F];      makeZero(h_Q,N_s*N_E*N_F);   

  scalar* h_Vtmp    = new scalar[N_s*N_E];          makeZero(h_Vtmp, N_s*N_E);
  scalar* h_dVinteg = new scalar[N_G*N_E];          makeZero(h_dVinteg, N_G*N_E);
  
  // copy from the fullMatrix to the pointer format (column major)
  copyMatrixToPointer(phi,h_phi);
  copyMatrixToPointer(phi_w,h_phi_w);
  copyMatrixToPointer(dphi,h_dphi);
  copyMatrixToPointer(dphi_w,h_dphi_w);
  for(int e = 0; e < N_E; e++){
    h_J[e] = J(e,0);
  }
  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
      for(int alpha = 0; alpha < D; alpha++){
  	for(int a = 0; a < D; a++){
  	  h_invJac[((e*N_G+g)*D+alpha)*D+a] = invJac(e*D+alpha,g*D+a);
  	}
      }
    }
  }
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      for(int i = 0; i < N_s; i++){
  	h_U[(e*N_F+fc)*N_s+i]= U(i,e*N_F+fc);
      }
    }
  }

  //
  // Initialize and allocate some stuff on the device
  //
  scalar* d_phi, *d_phi_w, *d_dphi, *d_dphi_w;
  scalar* d_J, *d_invJac;
  scalar* d_Minv;
  scalar* d_U, *d_Us, *d_Ustar, *d_DU, *d_UF;
  scalar* d_Uinteg, *d_dUinteg, *d_UintegF;
  scalar* d_s, *d_f, *d_q; 
  scalar* d_sJ, *d_fJ; 
  scalar* d_S, *d_F, *d_Q;

  if (!cpu){
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_phi,N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_phi_w,N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_dphi,D*N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_dphi_w,D*N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_J,N_E*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_invJac,N_G*D*N_E*D*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_Minv,N_s*N_s*N_E*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_U,N_s*N_E*N_F*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_Us,N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_Ustar,N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_DU,N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_UF,M_s*M_T*N_F*2*sizeof(scalar)));

    CUDA_SAFE_CALL(cudaMalloc((void**) &d_Uinteg,N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_dUinteg,D*N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_UintegF,M_G*M_T*N_F*2*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_s,N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_sJ,N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_S,N_s*N_E*N_F*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_f,D*N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_fJ,D*N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_F,N_s*N_E*N_F*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_q,M_G*M_T*N_F*2*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_Q,N_s*N_E*N_F*sizeof(scalar)));
  }

  //
  // Send the stuff to the device
  //
  if (!cpu){     
    CUDA_SAFE_CALL(cudaMemcpy(d_phi, h_phi, N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_phi_w, h_phi_w, N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_dphi, h_dphi, D*N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_dphi_w, h_dphi_w, D*N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_J, h_J, N_E*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_invJac, h_invJac, N_G*D*N_E*D*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_Minv, h_Minv, N_s*N_s*N_E*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_U, h_U, N_s*N_E*N_F*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemset(d_Q, (scalar)0.0, N_E*N_F*N_s*sizeof(scalar)))
  }

  
  //////////////////////////////////////////////////////////////////////////   
  //
  // Solve the problem on the CPU/GPU.
  //
  //////////////////////////////////////////////////////////////////////////   

  double T = 0;
  double Tstar = 0;
  scalar Dt = inputs.getTimeStep(); 
  int N_t = inputs.getNumberStep(); 
  int output_factor = inputs.getOutputFactor();
  int count = 1;
 
  // Runge-Kutta 4 coefficients
  scalar *beta = new scalar[4];
  beta[0] = 0.0; beta[1] = 0.5; beta[2] = 0.5; beta[3] = 1.0;
  scalar *gamma = new scalar[4];
  gamma[0] = 1.0/6.0; gamma[1] = 2.0/6.0; gamma[2] = 2.0/6.0; gamma[3] = 1.0/6.0;

  // Initialize the first DU evaluation
  // map U onto UF: requires Map, Ustar, UF and some integers for sizes, etc
  printf("=========== First DU evaluation ==============\n");
  if (cpu){
    if(multifluid)   Lcpu_mapToFace_multifluid(M_s, M_T, N_F, N_s, boundaryMap, h_U, h_UF);
    else if(passive) Lcpu_mapToFace_passive(M_s, M_T, N_F, N_s, boundaryMap, h_U, h_UF);
  }
  else if(!cpu){
    if(multifluid) Lgpu_mapToFace_multifluid(M_s, M_T, N_F, N_s, boundaryMap, d_U, d_UF);
    CUDA_SAFE_CALL(cudaThreadSynchronize());
  }
 
  if (debug){
    printf("Check mapping of U onto UF\n");
    if(!cpu){
      CUDA_SAFE_CALL(cudaMemcpy(h_U, d_U, N_F*N_s*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
      CUDA_SAFE_CALL(cudaMemcpy(h_UF, d_UF, M_T*N_F*2*sizeof(scalar), cudaMemcpyDeviceToHost));
    }
    for(int e = 0; e < N_E; e++){
      int fc = 2;
      printf("%f ?= %f    and    ", h_U[(e*N_F+fc)*N_s+0], h_UF[(e*N_F+fc)*2+1]);
      printf("%f ?= %f\n", h_U[(e*N_F+fc)*N_s+1], h_UF[((e+1)*N_F+fc)*2+0]);
    }
  }

  // Get the velocity field (to later find the derivative)
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      h_Vtmp[e*N_s+i] = h_U[(e*N_F+1)*N_s+i]/h_U[(e*N_F+0)*N_s+i]; // u = (rho u)/rho
    }
  }
  
  // collocationU: requires phi, dphi, Ustar, Uinteg, dUinteg and some sizes
  if (cpu){
    if (blas == 1){
      blasGemm("N","N", N_G  , N_E*N_F, N_s, 1, h_phi,  N_G  , h_U, N_s, 0.0, h_Uinteg , N_G);
      blasGemm("N","N", N_G*D, N_E*N_F, N_s, 1, h_dphi, N_G*D, h_U, N_s, 0.0, h_dUinteg, N_G*D);
      blasGemm("N","N", N_G, N_E, N_s, 1, h_dphi, N_G, h_Vtmp, N_s, 0.0, h_dVinteg, N_G);}
    else Lcpu_collocationU(D, N_G, N_s, N_E, N_F, h_Uinteg, h_dUinteg, h_phi, h_dphi, h_U);
  }
  else if(!cpu){
    if (blas == 1){
      cublasGemm('N','N', N_G  , N_E*N_F, N_s, 1, d_phi , N_G  , d_U, N_s, 0.0, d_Uinteg , N_G);
      cublasGemm('N','N', N_G*D, N_E*N_F, N_s, 1, d_dphi, N_G*D, d_U, N_s, 0.0, d_dUinteg, N_G*D);}
    else Lgpu_collocationU(D, N_G, N_s, N_E, N_F, d_Uinteg, d_dUinteg, d_phi, d_dphi, d_U);
    CUDA_SAFE_CALL(cudaThreadSynchronize());
  }
  
  // collocationUF: requires psi, UF, UintegF and some sizes
  if (cpu){
    for(int k = 0; k < 2*N_F*M_T; k++){ h_UintegF[k] = h_UF[k];}
  }
  else if(!cpu){
    cublasCopy(2*N_F*M_T, d_UF, 1, d_UintegF, 1);
    CUDA_SAFE_CALL(cudaThreadSynchronize());
  }
  
  if (debug){
    printf("Check collocation of U\n");
    if(!cpu) CUDA_SAFE_CALL(cudaMemcpy(h_Uinteg, d_Uinteg, N_F*N_G*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
    printf("===== test Uinteg  =====\n");
    for(int g = 0; g < N_G; g++){
      for(int e = 0; e < N_E; e++){
  	for(int fc = 0; fc < N_F; fc++){
  	  printf("%f ",h_Uinteg[(e*N_F+fc)*N_G+g]);
  	}
      }
      printf("\n");
    }
    printf("=== end test Uinteg ==\n");
    if(!cpu) CUDA_SAFE_CALL(cudaMemcpy(h_UintegF, d_UintegF, 2*N_F*M_G*M_T*sizeof(scalar), cudaMemcpyDeviceToHost));
    printf("===== test UintegF GPU eta =====\n");
    for(int g = 0; g < M_G; g++){
      for(int t = 0; t < M_T; t++){
    	printf("%4.3f %4.3f | ", h_UintegF[((t*N_F+0)*2+0)*M_G+g] , h_UintegF[((t*N_F+0)*2+1)*M_G+g]);
      }
      printf("\n");
    }
    printf("=== end test UintegF GPU eta ==\n");
  }

  // evaluate_sf: requires Uinteg, (dUintegR), H0, G0, s,f 
  if (cpu){
    if(multifluid) Lcpu_evaluate_sf_multifluid(D, N_G, N_E, N_F, model, h_s, h_f, h_Uinteg, h_dUinteg, h_invJac);
    if(passive)    Lcpu_evaluate_sf_passive(D, N_G, N_E, N_F, gamma0, h_s, h_f, h_Uinteg, h_dUinteg, h_invJac);
    //if(multifluid) Lcpu_evaluate_sf_multifluid(D, N_G, N_E, N_F, model, h_s, h_f, h_Uinteg, h_dVinteg, h_invJac);
  }
  else if(!cpu){
    if(multifluid) Lgpu_evaluate_sf_multifluid(D, N_G, N_E, N_F, model, d_s, d_f, d_Uinteg, d_dUinteg, d_invJac);
    CUDA_SAFE_CALL(cudaThreadSynchronize());
  }
    
  if (debug){
    if(!cpu) CUDA_SAFE_CALL(cudaMemcpy(h_f, d_f, D*N_F*N_G*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
    printf("===== test f GPU eta =====\n");
    for(int e = 0; e < N_E; e++){
      for(int g = 0; g < N_G; g++){
  	printf("element %i, eta: f_x=%f \n", e, h_f[((e*N_F+0)*N_G+g)*D+0]);
  	printf("element %i,  ux: f_x=%f \n", e, h_f[((e*N_F+1)*N_G+g)*D+0]);
  	printf("element %i,  uy: f_x=%f \n", e, h_f[((e*N_F+2)*N_G+g)*D+0]);
      }
    }
    printf("=== end test f GPU eta ==\n");
  }

  // evaluate_q: requires UintegF, normals, q, H0, G0
  if (cpu){
    if(multifluid) Lcpu_evaluate_q_multifluid(M_G, M_T, N_F, flux, model, h_q, h_UintegF);
    if(passive)    Lcpu_evaluate_q_passive(M_G, M_T, N_F, flux, gamma0, h_q, h_UintegF);
  }
  else if (!cpu){
    if(multifluid) Lgpu_evaluate_q_multifluid(M_G, M_T, N_F, model, d_q, d_UintegF);
    CUDA_SAFE_CALL(cudaThreadSynchronize());
  }
  
  // redistribute_sf: requires J, invJac, s, f, phi_w, dphi_w, sJ, fJ, S, F
  if (cpu){
    Lcpu_redistribute_sf(D, N_G, N_E, N_F, h_sJ, h_fJ, h_s, h_f, h_J, h_invJac);
  }
  else if (!cpu){
    Lgpu_redistribute_sf(D, N_G, N_E, N_F, d_sJ, d_fJ, d_s, d_f, d_J, d_invJac);
    CUDA_SAFE_CALL(cudaThreadSynchronize());
  }

  if (debug){
    if(!cpu){
      CUDA_SAFE_CALL(cudaMemcpy(h_sJ, d_sJ, N_F*N_G*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
      CUDA_SAFE_CALL(cudaMemcpy(h_fJ, d_fJ, D*N_F*N_G*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
    }
    // test sJ
    printf("test sJ GPU\n");
    for (int e = 0; e < N_E; e++){
      for (int g = 0; g < N_G; g++){
  	for (int fc = 0; fc < N_F; fc++){
  	  printf("%7.4f", h_sJ[(e*N_F+fc)*N_G+g]);
  	}
      }
      printf("\n");
    }
    printf("end test sJ GPU\n");
    // test fJ
    printf("test fJ GPU\n");
    for (int e = 0; e < N_E; e++){
      for (int g = 0; g < N_G; g++){
  	for (int fc = 0; fc < N_F; fc++){
	  printf("%7.4f", h_fJ[((e*N_F+fc)*N_G+g)*D+0]);
    	}
      }
      printf("\n");
    }
    printf("end test fJ GPU\n");
  }

  // matrix-matrix for sf
  if (cpu){
    if (blas==1) {
      blasGemm("T","N", N_s, N_E*N_F, N_G  , 1, h_phi_w , N_G  , h_sJ, N_G  , 0.0, h_S, N_s);
      blasGemm("T","N", N_s, N_E*N_F, N_G*D, 1, h_dphi_w, N_G*D, h_fJ, N_G*D, 0.0, h_F, N_s);
    }
    else Lcpu_gemm_sf(D, N_G, N_s, N_E, N_F, h_S, h_F, h_sJ, h_fJ, h_phi_w, h_dphi_w);
  }
  else if (!cpu){
    if (blas==1) {
      cublasGemm('T','N', N_s, N_E*N_F, N_G  , 1, d_phi_w , N_G  , d_sJ, N_G  , 0.0, d_S, N_s);
      cublasGemm('T','N', N_s, N_E*N_F, N_G*D, 1, d_dphi_w, N_G*D, d_fJ, N_G*D, 0.0, d_F, N_s);
    }
    else Lgpu_gemm_sf(D, N_G, N_s, N_E, N_F, d_S, d_F, d_sJ, d_fJ, d_phi_w, d_dphi_w);
  }
  

  // map_q: requires map, Qtcj, Q (might want to do this in the previous step)
  if (cpu){
    Lcpu_mapToElement(N_s, N_E, N_F, h_Q, h_q);
  }
  else if (!cpu){
    Lgpu_mapToElement(N_s, N_E, N_F, d_Q, d_q);
    CUDA_SAFE_CALL(cudaThreadSynchronize());
  }

  if (debug){
    if(!cpu){
      CUDA_SAFE_CALL(cudaMemcpy(h_S, d_S, N_F*N_s*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
      CUDA_SAFE_CALL(cudaMemcpy(h_F, d_F, N_F*N_s*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
      CUDA_SAFE_CALL(cudaMemcpy(h_Q, d_Q, N_F*N_s*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
    }
    // Test S, F, Q
    printf("============ test GPU SFQ ===============\n");
    for(int i = 0; i < N_s; i++){
      for(int e = 0; e < N_E; e++){
  	printf("s:%7.4f f:%7.4f q:%7.4f", h_S[(e*N_F+0)*N_s+i], h_F[(e*N_F+0)*N_s+i], h_Q[(e*N_F+0)*N_s+i]);
      }
      printf("\n");
    }
    printf("============ end test GPU SFQ ===========\n");
    printf("============ test GPU SFQx ==============\n");
    for(int i = 0; i < N_s; i++){
      for(int e = 0; e < N_E; e++){
  	printf("s:%7.4f f:%7.4f q:%7.4f", h_S[(e*N_F+1)*N_s+i], h_F[(e*N_F+1)*N_s+i], h_Q[(e*N_F+1)*N_s+i]);
      }
      printf("\n");
    }
    printf("============ end test GPU SFQx ===========\n");
    printf("============ test GPU SFQy ==============\n");
    for(int i = 0; i < N_s; i++){
      for(int e = 0; e < N_E; e++){
  	printf("s:%7.4f f:%7.4f q:%7.4f", h_S[(e*N_F+2)*N_s+i], h_F[(e*N_F+2)*N_s+i], h_Q[(e*N_F+2)*N_s+i]);
      }
      printf("\n");
    }
    printf("============ end test GPU SFQy ===========\n");
  }

  // solve: requires Q, F, S, Dt, Minv, DU
  if (cpu){
    Lcpu_solve(N_s, N_E, N_F, h_DU, h_S, h_F, h_Q, h_Minv, Dt);
  }
  else if (!cpu){
    Lgpu_solve(N_s, N_E, N_F, d_DU, d_S, d_F, d_Q, d_Minv, Dt);
    CUDA_SAFE_CALL(cudaThreadSynchronize());
  }

  if (debug){
    if(!cpu) CUDA_SAFE_CALL(cudaMemcpy(h_DU, d_DU, N_F*N_s*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
    // Test DU
    printf("============ test GPU DU ===============\n");
    for(int i = 0; i < N_s; i++){
      for(int e = 0; e < N_E; e++){
  	printf("%7.4f ", h_DU[(e*N_F+0)*N_s+i]);
      }
      printf("\n");
    }
    printf("============ end test GPU DU ===========\n");
    printf("============ test GPU DU x ===============\n");
    for(int i = 0; i < N_s; i++){
      for(int e = 0; e < N_E; e++){
  	printf("%7.4f ", h_DU[(e*N_F+1)*N_s+i]);
      }
      printf("\n");
    }
    printf("============ end test GPU DU x ===========\n");    
    printf("============ test GPU DU y ===============\n");
    for(int i = 0; i < N_s; i++){
      for(int e = 0; e < N_E; e++){
  	printf("%7.4f ", h_DU[(e*N_F+2)*N_s+i]);
      }
      printf("\n");
    }
    printf("============ end test GPU DU y ===========\n");
  }

  
  // if 0-order average the solution in the cells
  if (order0){
    if (cpu){
      Lcpu_average_cell_p0(N_s, N_E, N_F, h_DU);
    }
    else if (!cpu){
      Lgpu_average_cell_p0(N_s, N_E, N_F, d_DU);
      CUDA_SAFE_CALL(cudaThreadSynchronize());
    }
  }

  printf("==== Now RK 4 steps =====\n");
  
  // Time  integration  
  for (int n = 1; n <= N_t; n++){

    //
    // RK4
    //
    CUDA_SAFE_CALL(cudaThreadSynchronize());
    if (cpu){
      if (blas==1) {blasCopy(N_F*N_s*N_E, h_U, 1, h_Us, 1);}    
      else Lcpu_equal(N_s, N_E, N_F, h_Us, h_U);// make Us = U;
    }
    else if (!cpu){
      if (blas==1) {cublasCopy(N_F*N_s*N_E, d_U, 1, d_Us, 1);}    
      else Lgpu_equal(N_s, N_E, N_F, d_Us, d_U);// make Us = U;
      CUDA_SAFE_CALL(cudaThreadSynchronize());
    }

    for(int k = 0; k < 4; k++){
      if (cpu){
  	if (blas==1) {blasCopy(N_F*N_s*N_E, h_Us, 1, h_Ustar, 1);}    
  	else Lcpu_equal(N_s, N_E, N_F, h_Ustar, h_Us); // make Ustar = Us;
      }
      else if (!cpu){
  	if (blas==1) {cublasCopy(N_F*N_s*N_E, d_Us, 1, d_Ustar, 1);}    
  	else Lgpu_equal(N_s, N_E, N_F, d_Ustar, d_Us); // make Ustar = Us;
  	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }
      if (cpu){
  	if (blas==1) {blasAxpy(N_s*N_F*N_E, beta[k], h_DU, 1, h_Ustar, 1);}      
  	else Lcpu_add(N_s, N_E, N_F, h_Ustar, h_DU, beta[k]);// do Ustar.add(DU,beta[k]);
      }
      else if (!cpu){
  	if (blas==1) {cublasAxpy(N_s*N_F*N_E, beta[k], d_DU, 1, d_Ustar, 1);}      
  	else Lgpu_add(N_s, N_E, N_F, d_Ustar, d_DU, beta[k]);// do Ustar.add(DU,beta[k]);
  	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }
      Tstar = T + beta[k]*Dt;

      // map U onto UF: requires Map, Ustar, UF and some integers for sizes, etc
      if (cpu){
	if(multifluid)   Lcpu_mapToFace_multifluid(M_s, M_T, N_F, N_s, boundaryMap, h_Ustar, h_UF);
	else if(passive) Lcpu_mapToFace_passive(M_s, M_T, N_F, N_s, boundaryMap, h_Ustar, h_UF);
      }
      else if(!cpu){
      	if(multifluid) Lgpu_mapToFace_multifluid(M_s, M_T, N_F, N_s, boundaryMap, d_Ustar, d_UF);
   	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }

      // Get the velocity field (to later find the derivative)
      for(int e = 0; e < N_E; e++){
	for(int i = 0; i < N_s; i++){
	  h_Vtmp[e*N_s+i] = h_Ustar[(e*N_F+1)*N_s+i]/h_Ustar[(e*N_F+0)*N_s+i]; // u = (rho u)/rho
	}
      }

      // collocationU: requires phi, dphi, Ustar, Uinteg, dUinteg and some sizes
      if (cpu){
  	if (blas==1) {
  	  blasGemm("N","N", N_G  , N_E*N_F, N_s, 1, h_phi,  N_G  , h_Ustar, N_s, 0.0, h_Uinteg , N_G);
  	  blasGemm("N","N", N_G*D, N_E*N_F, N_s, 1, h_dphi, N_G*D, h_Ustar, N_s, 0.0, h_dUinteg, N_G*D);
	  blasGemm("N","N", N_G, N_E, N_s, 1, h_dphi, N_G, h_Vtmp, N_s, 0.0, h_dVinteg, N_G);}
  	else Lcpu_collocationU(D, N_G, N_s, N_E, N_F, h_Uinteg, h_dUinteg, h_phi, h_dphi, h_Ustar);
      }
      else if(!cpu){
  	if (blas==1) {
  	  cublasGemm('N','N', N_G  , N_E*N_F, N_s, 1, d_phi , N_G  , d_Ustar, N_s, 0.0, d_Uinteg , N_G);
  	  cublasGemm('N','N', N_G*D, N_E*N_F, N_s, 1, d_dphi, N_G*D, d_Ustar, N_s, 0.0, d_dUinteg, N_G*D);}
  	else Lgpu_collocationU(D, N_G, N_s, N_E, N_F, d_Uinteg, d_dUinteg, d_phi, d_dphi, d_Ustar);
  	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }

      // collocationUF: requires psi, UF, UintegF and some sizes
      if (cpu){
	for(int k = 0; k < 2*N_F*M_T; k++){ h_UintegF[k] = h_UF[k];}
      }
      else if(!cpu){
	cublasCopy(2*N_F*M_T, d_UF, 1, d_UintegF, 1);
	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }
     
      // evaluate_sf: requires Uinteg, (dUintegR), H0, G0, s,f 
      if (cpu){
	if(multifluid) Lcpu_evaluate_sf_multifluid(D, N_G, N_E, N_F, model, h_s, h_f, h_Uinteg, h_dUinteg, h_invJac);
	if(passive)    Lcpu_evaluate_sf_passive(D, N_G, N_E, N_F, gamma0, h_s, h_f, h_Uinteg, h_dUinteg, h_invJac);
	//if(multifluid) Lcpu_evaluate_sf_multifluid(D, N_G, N_E, N_F, model, h_s, h_f, h_Uinteg, h_dVinteg, h_invJac);
      }
      else if(!cpu){
	if(multifluid) Lgpu_evaluate_sf_multifluid(D, N_G, N_E, N_F, model, d_s, d_f, d_Uinteg, d_dUinteg, d_invJac);
  	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }

      // evaluate_q: requires UintegF, normals, q, H0, G0
      if (cpu){
	if(multifluid) Lcpu_evaluate_q_multifluid(M_G, M_T, N_F, flux, model, h_q, h_UintegF);
	if(passive)    Lcpu_evaluate_q_passive(M_G, M_T, N_F, flux, gamma0, h_q, h_UintegF);
      }
      else if (!cpu){
	if(multifluid) Lgpu_evaluate_q_multifluid(M_G, M_T, N_F, model, d_q, d_UintegF);
  	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }

      // redistribute_sf: requires J, invJac, s, f, phi_w, dphi_w, sJ, fJ, S, F
      if (cpu){
   	Lcpu_redistribute_sf(D, N_G, N_E, N_F, h_sJ, h_fJ, h_s, h_f, h_J, h_invJac);
      }
      else if (!cpu){
  	Lgpu_redistribute_sf(D, N_G, N_E, N_F, d_sJ, d_fJ, d_s, d_f, d_J, d_invJac);
  	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }

      // matrix-matrix multiply for sf
      if (cpu){
	if (blas==1)  {
	  blasGemm("T","N", N_s, N_E*N_F, N_G  , 1, h_phi_w , N_G  , h_sJ, N_G  , 0.0, h_S, N_s);
	  blasGemm("T","N", N_s, N_E*N_F, N_G*D, 1, h_dphi_w, N_G*D, h_fJ, N_G*D, 0.0, h_F, N_s);}
	else Lcpu_gemm_sf(D, N_G, N_s, N_E, N_F, h_S, h_F, h_sJ, h_fJ, h_phi_w, h_dphi_w);
      }
      else if (!cpu){
      	if (blas==1)  {
      	  cublasGemm('T','N', N_s, N_E*N_F, N_G  , 1, d_phi_w , N_G  , d_sJ, N_G  , 0.0, d_S, N_s);
      	  cublasGemm('T','N', N_s, N_E*N_F, N_G*D, 1, d_dphi_w, N_G*D, d_fJ, N_G*D, 0.0, d_F, N_s);}
      	else Lgpu_gemm_sf(D, N_G, N_s, N_E, N_F, d_S, d_F, d_sJ, d_fJ, d_phi_w, d_dphi_w);
      	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }

      
      // map_q: requires map, Qtcj, Q (might want to do this in the previous step)
      if (cpu){
	Lcpu_mapToElement(N_s, N_E, N_F, h_Q, h_q);
      }
      else if (!cpu){
	Lgpu_mapToElement(N_s, N_E, N_F, d_Q, d_q);
      	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }

      // solve: requires Q, F, S, Dt, Minv, DU 
      if (cpu){
	Lcpu_solve(N_s, N_E, N_F, h_DU, h_S, h_F, h_Q, h_Minv, Dt);
      }
      else if (!cpu){
      	Lgpu_solve(N_s, N_E, N_F, d_DU, d_S, d_F, d_Q, d_Minv, Dt);
      	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }


      // if 0-order average the solution in the cells
      if (order0){
	if (cpu){
	  Lcpu_average_cell_p0(N_s, N_E, N_F, h_DU);
	}
	else if (!cpu){
	  Lgpu_average_cell_p0(N_s, N_E, N_F, d_DU);
	  CUDA_SAFE_CALL(cudaThreadSynchronize());
	}
      }
      
      if (cpu){
	if (blas==1) {blasAxpy(N_s*N_F*N_E, gamma[k], h_DU, 1, h_U, 1);}      
	else Lcpu_add(N_s, N_E, N_F, h_U, h_DU, gamma[k]); // do U.add(DU,gamma[k])
      }
      else if (!cpu){
      	if (blas==1) {cublasAxpy(N_s*N_F*N_E, gamma[k], d_DU, 1, d_U, 1);}      
      	else Lgpu_add(N_s, N_E, N_F, d_U, d_DU, gamma[k]); // do U.add(DU,gamma[k]
      	CUDA_SAFE_CALL(cudaThreadSynchronize());
      }
    } // end RK4 loop


    T = T + Dt;


    //
    // Get the solution on the CPU so that we can 
    // output it to a file 
    // 

    if(n % (N_t/output_factor) == 0){

      // Get the solution to the CPU
      if (!cpu){
      	CUDA_SAFE_CALL(cudaMemcpy(h_U, d_U, N_s*N_F*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
      }
      
      printf("Solution written to output file at step %i and time %f.\n",n,n*Dt);
      if(multifluid)print_dg_multifluid(N_s, N_E, N_F, model, h_U, m, msh_lin, count, n*Dt, 1,-1);
      if(passive)   print_dg_passive(N_s, N_E, N_F, gamma0, h_U, m, msh_lin, count, n*Dt, 1,-1);
      count++;
    }
        
  } // end time integration

    

  //////////////////////////////////////////////////////////////////////////   
  //
  // Free stuff on the device
  //
  //////////////////////////////////////////////////////////////////////////   
  status = cublasShutdown();
  if (!cpu){
    CUDA_SAFE_CALL(cudaFree(d_phi));
    CUDA_SAFE_CALL(cudaFree(d_phi_w));
    CUDA_SAFE_CALL(cudaFree(d_dphi));
    CUDA_SAFE_CALL(cudaFree(d_dphi_w));
    CUDA_SAFE_CALL(cudaFree(d_J));
    CUDA_SAFE_CALL(cudaFree(d_invJac));
    CUDA_SAFE_CALL(cudaFree(d_Minv));
    CUDA_SAFE_CALL(cudaFree(d_U));
    CUDA_SAFE_CALL(cudaFree(d_Us));
    CUDA_SAFE_CALL(cudaFree(d_Ustar));
    CUDA_SAFE_CALL(cudaFree(d_DU));
    CUDA_SAFE_CALL(cudaFree(d_UF));
    CUDA_SAFE_CALL(cudaFree(d_Uinteg));
    CUDA_SAFE_CALL(cudaFree(d_dUinteg));
    CUDA_SAFE_CALL(cudaFree(d_UintegF));
    CUDA_SAFE_CALL(cudaFree(d_s));
    CUDA_SAFE_CALL(cudaFree(d_sJ));
    CUDA_SAFE_CALL(cudaFree(d_S));
    CUDA_SAFE_CALL(cudaFree(d_f));
    CUDA_SAFE_CALL(cudaFree(d_fJ));
    CUDA_SAFE_CALL(cudaFree(d_F));
    CUDA_SAFE_CALL(cudaFree(d_q));
    CUDA_SAFE_CALL(cudaFree(d_Q));
  }
  

  //////////////////////////////////////////////////////////////////////////   
  //
  // Output some stuff in a file to read by post-proc
  //
  //////////////////////////////////////////////////////////////////////////
  std::string post = "post.dat"; 
  FILE *f = fopen(post.c_str(),"w");
  fprintf(f,"%12.5E\n", XYZNodes(1, 0*D+0)-XYZNodes(0, 0*D+0));


  
  //////////////////////////////////////////////////////////////////////////   
  //
  // Free stuff on the host
  //
  //////////////////////////////////////////////////////////////////////////   

  delete[] h_Minv;
  delete[] h_phi;
  delete[] h_phi_w;
  delete[] h_dphi;
  delete[] h_dphi_w;
  delete[] h_J;
  delete[] h_invJac;
  delete[] h_U;
  delete[] h_Us;
  delete[] h_Ustar;
  delete[] h_DU;
  delete[] h_UF;
  delete[] h_Uinteg;
  delete[] h_dUinteg;
  delete[] h_UintegF;
  delete[] h_s;
  delete[] h_sJ;
  delete[] h_S;
  delete[] h_f;
  delete[] h_fJ;
  delete[] h_F;
  delete[] h_q;
  delete[] h_Q;

  delete[] h_Vtmp;
  delete[] h_dVinteg;
  
  delete[] beta;
  delete[] gamma;

}// end main



//
// Function definitions
//
void copyMatrixToPointer(fullMatrix<scalar> &A, scalar* h_A){
  
  // Column major sorting
  for(int j = 0; j < A.size2(); j++){
    for(int i = 0; i < A.size1(); i++){
      h_A[j*A.size1()+i] = A(i,j);
    }
  }
}

void copyMatrixToPointer(const fullMatrix<scalar> &A, scalar* h_A){
  
  // Column major sorting
  for(int j = 0; j < A.size2(); j++){
    for(int i = 0; i < A.size1(); i++){
      h_A[j*A.size1()+i] = A(i,j);
    }
  }
}
void average_cell_p0(const int N_s, const int N_E, const int N_F, fullMatrix<scalar> &U){
  fullMatrix<scalar> average(N_s,N_s);
  average.setAll((scalar)1.0/N_s);
  fullMatrix<scalar> tmp(N_s,N_E*N_F);
  tmp.gemm(average,U);
  U=tmp; 
}
