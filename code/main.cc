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
#include "polynomialsJacobi.h"
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
#include <limiter_fct.h>

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

void vandermonde1d(const int order, const fullMatrix<scalar> r, fullMatrix<scalar> &V1D);
void monovandermonde1d(const int order, const fullMatrix<scalar> r, fullMatrix<scalar> &V1D);

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
  if      (inputs.getFlux()=="llf") {flux = 0; printf("Using LLF\n");}
  else if (inputs.getFlux()=="ncf") {flux = 1; printf("Using NCF\n");}
  else if (inputs.getFlux()=="roe") {flux = 2; printf("Using ROE\n");}
  else{ printf("Invalid flux setup. Correct the deck.\n");}
  
  // Get the mesh
  std::string fileName = inputs.getMeshfile();
  
  // Setup the problem type
  bool multifluid = false;
  bool passive = false;
  if      (inputs.getProblem() == "multifluid")   multifluid = true;
  else if (inputs.getProblem() == "passive")      passive = true;
  else{ printf("Invalid problem setup. Correct the deck.\n");}

  // Setup the model
  int model = 0;
  if      (inputs.getModel() == "gammamod")   model = 0;
  else if (inputs.getModel() == "invgamma")   model = 1;
  else if (passive) model = 0;
  else{ printf("Invalid model setup. Correct the deck.\n");}

  // Setup the limiting
  int limiter = 0;
  if      (inputs.getLimiter() == "hsl")   {limiter = 1; printf("Using HS limiting\n");}
  else if (inputs.getLimiter() == "hrl")   {limiter = 2; printf("Using HR limiting\n");}
  else{limiter = 0; printf("No limiting\n");}
  
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

  //////////////////////////////////////////////////////////////////////////   
  //
  // Monomial to Lagrange basis transforms
  //
  //////////////////////////////////////////////////////////////////////////
  fullMatrix<scalar> points1D(N_G,1); for(int g=0; g<N_G; g++) {points1D(g,0) = points(g,0);}
  fullMatrix<scalar> phiinv(N_s,N_G);
  phi.invert(phiinv);
  fullMatrix<scalar> monoV1D;
  fullMatrix<scalar> monoV1Dinv;
  monovandermonde1d(order, points1D, monoV1D);
  monoV1D.invert(monoV1Dinv);

  // for(int i = 0; i<monoV1D.size1();i++){
  //   for(int j = 0; j<monoV1D.size2();j++){
  //     printf("%f\t",monoV1D(i,j));
  //   }
  //   printf("\n");
  // }

  // Go from lagrange to monomial basis 
  fullMatrix<scalar> Lag2Mono(N_s,N_s);
  fullMatrix<scalar> Mono2Lag(N_s,N_s);
  Lag2Mono.gemm(monoV1Dinv, phi);   // Calculate the complete nodal to modal transform = V1Dinv*phiGL
  Lag2Mono.invert(Mono2Lag);
  //Mono2Lag.gemm(phiinv,monoV1D);
  

  //////////////////////////////////////////////////////////////////////////   
  //
  // Modal basis functions and change of basis
  //
  //////////////////////////////////////////////////////////////////////////

  // Get the Legendre-Gauss-Lobatto points
  fullMatrix<scalar> pointsGL;
  JacobiGL(0,0,order,pointsGL);
  int N_GL = pointsGL.size1();
  
  // Evaluate the Lagrange polynomials at the GL points (needed for collocation)
  fullMatrix<scalar> phiGL   (N_GL,N_s);
  fullMatrix<scalar> phiGLinv(N_s,N_GL); 
  fullMatrix<double> phiDGL  (N_GL,N_s); 
  basis->f (pointsGL, phiDGL);
  for(int gl = 0; gl < N_GL; gl++){
    for(int i = 0; i < N_s; i++){
      phiGL(gl,i) = (scalar)phiDGL(gl,i);
    }	  
  }   
  phiGL.invert(phiGLinv);
  
  // Vandermonde matrix to go from modal to nodal representation
  // u_nodal = V1D    u_modal
  // u_modal = V1Dinv u_nodal
  fullMatrix<scalar> V1D;
  fullMatrix<scalar> V1Dinv;
  vandermonde1d(order, pointsGL, V1D);
  V1D.invert(V1Dinv);
  
  // Go from nodal to modal
  fullMatrix<scalar> Nod2Mod(N_s,N_s);
  fullMatrix<scalar> Mod2Nod(N_s,N_s);
  Nod2Mod.gemm(V1Dinv, phiGL);   // Calculate the complete nodal to modal transform = V1Dinv*phiGL
  Nod2Mod.invert(Mod2Nod);
  
  
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


  // // Setup an dummy initial condition to validate the limiting
  // if(order==1){
  //   for(int fc=0; fc<N_F; fc++){
  //     U(0,0*N_F+fc) = 1;
  //     U(1,0*N_F+fc) = 2;
  //     U(0,1*N_F+fc) = 3;
  //     U(1,1*N_F+fc) = 1.5;
  //     U(0,2*N_F+fc) = 1;
  //     U(1,2*N_F+fc) = 0.5;
  //   }
  // }
  // else if (order==2){
  //   printf("hello\n");
  //   for(int fc=0; fc<N_F; fc++){
  //     U(0,0*N_F+fc) = 1;
  //     U(1,0*N_F+fc) = 1.5;
  //     U(2,0*N_F+fc) = 0.75;
  //     U(0,1*N_F+fc) = 1.75;
  //     U(1,1*N_F+fc) = 0.5;
  //     U(2,1*N_F+fc) = 1.25;
  //     U(0,2*N_F+fc) = 1;
  //     U(1,2*N_F+fc) = 0.75;
  //     U(2,2*N_F+fc) = 2;
  //   }
  // }

  
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
#ifdef USE_GPU
  // Choose the device
  CUDA_SAFE_CALL(cudaSetDevice(0));

  // Use cublas or not
  cublasStatus status;
  status = cublasInit();
#endif
  
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
  scalar* h_phiGL   = new scalar[N_GL*N_s];         makeZero(h_phiGL,N_GL*N_s);
  scalar* h_phiGLinv= new scalar[N_GL*N_s];         makeZero(h_phiGLinv,N_GL*N_s);
  scalar* h_phi_w   = new scalar[N_G*N_s];          makeZero(h_phi_w,N_G*N_s);          
  scalar* h_dphi    = new scalar[D*N_G*N_s];	    makeZero(h_dphi,D*N_G*N_s);	 
  scalar* h_dphi_w  = new scalar[D*N_G*N_s];	    makeZero(h_dphi_w,D*N_G*N_s);
  scalar* h_V1D     = new scalar[N_GL*N_s];         makeZero(h_V1D,N_GL*N_s);
  scalar* h_V1Dinv  = new scalar[N_GL*N_s];         makeZero(h_V1Dinv,N_s*N_GL);
  scalar* h_monoV1D     = new scalar[N_GL*N_s];         makeZero(h_monoV1D,N_GL*N_s);
  scalar* h_monoV1Dinv  = new scalar[N_GL*N_s];         makeZero(h_monoV1Dinv,N_s*N_GL);
  scalar* h_weight  = new scalar[N_G];              makeZero(h_weight,N_G);
  scalar* h_J       = new scalar[N_E];              makeZero(h_J,N_E);                                 // not same as J!!
  scalar* h_invJac  = new scalar[N_G*D*N_E*D];      makeZero(h_invJac,N_G*D*N_E*D);                    // not same as invJac!!
  scalar* h_Us      = new scalar[N_s*N_E*N_F];	    makeZero(h_Us,N_s*N_E*N_F);	 
  scalar* h_Ustar   = new scalar[N_s*N_E*N_F];	    makeZero(h_Ustar,N_s*N_E*N_F);
  scalar* h_DU      = new scalar[N_s*N_E*N_F];	    makeZero(h_DU,N_s*N_E*N_F);	 
  scalar* h_U       = new scalar[N_s*N_E*N_F];	    makeZero(h_U,N_s*N_E*N_F);
  scalar* h_UMod    = new scalar[N_s*N_E*N_F];      makeZero(h_UMod,N_s*N_E*N_F);
  scalar* h_UModNew = new scalar[N_s*N_E*N_F];      makeZero(h_UModNew,N_s*N_E*N_F);
  scalar* h_A       = new scalar[N_s*N_E*N_F];      makeZero(h_A,N_s*N_E*N_F);
  scalar* h_Alim    = new scalar[N_s*N_E*N_F];      makeZero(h_Alim,N_s*N_E*N_F);
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

  scalar* h_Nod2Mod = new scalar[N_s*N_s];          makeZero(h_Nod2Mod, N_s*N_s);
  scalar* h_Mod2Nod = new scalar[N_s*N_s];          makeZero(h_Mod2Nod, N_s*N_s);

  scalar* h_Lag2Mono = new scalar[N_s*N_s];         makeZero(h_Lag2Mono, N_s*N_s);
  scalar* h_Mono2Lag = new scalar[N_s*N_s];         makeZero(h_Mono2Lag, N_s*N_s);
  
  scalar* h_Vtmp    = new scalar[N_s*N_E];          makeZero(h_Vtmp, N_s*N_E);
  scalar* h_dVinteg = new scalar[N_G*N_E];          makeZero(h_dVinteg, N_G*N_E);
  
  // copy from the fullMatrix to the pointer format (column major)
  copyMatrixToPointer(phi,h_phi);
  copyMatrixToPointer(phiGL,h_phiGL);
  copyMatrixToPointer(phiGLinv,h_phiGLinv);
  copyMatrixToPointer(V1D,h_V1D);
  copyMatrixToPointer(V1Dinv,h_V1Dinv);
  copyMatrixToPointer(monoV1D,h_monoV1D);
  copyMatrixToPointer(monoV1Dinv,h_monoV1Dinv);
  for(int g=0; g<N_G; g++) h_weight[g] = weight(g,0);
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
  copyMatrixToPointer(Nod2Mod,h_Nod2Mod);
  copyMatrixToPointer(Mod2Nod,h_Mod2Nod);
  copyMatrixToPointer(Lag2Mono,h_Lag2Mono);
  copyMatrixToPointer(Mono2Lag,h_Mono2Lag);
  
  
  //
  // Initialize and allocate some stuff on the device
  //
#ifdef USE_GPU
  scalar* d_phi, *d_phi_w, *d_dphi, *d_dphi_w;
  scalar* d_J, *d_invJac;
  scalar* d_Minv;
  scalar* d_U, *d_Us, *d_Ustar, *d_DU, *d_UF;
  scalar* d_Uinteg, *d_dUinteg, *d_UintegF;
  scalar* d_s, *d_f, *d_q; 
  scalar* d_sJ, *d_fJ; 
  scalar* d_S, *d_F, *d_Q;

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
#endif

  //
  // Send the stuff to the device
  //
#ifdef USE_GPU
  CUDA_SAFE_CALL(cudaMemcpy(d_phi, h_phi, N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_phi_w, h_phi_w, N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_dphi, h_dphi, D*N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_dphi_w, h_dphi_w, D*N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_J, h_J, N_E*sizeof(scalar), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_invJac, h_invJac, N_G*D*N_E*D*sizeof(scalar), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_Minv, h_Minv, N_s*N_s*N_E*sizeof(scalar), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_U, h_U, N_s*N_E*N_F*sizeof(scalar), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemset(d_Q, (scalar)0.0, N_E*N_F*N_s*sizeof(scalar)))
#endif

  
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

      // Limit the solution if you so want to do so
      if(k>0){
      	if(cpu){
      	  if (limiter==2){ //HR limiting
	    // Go from conservative to primitive space
	    //Lcpu_Cons2Prim(N_s, N_E, N_F, h_Ustar, multifluid, passive, model, gamma0);
      	    // Go from lagrange to monomial representation
      	    blasGemm("N","N", N_s, N_E*N_F, N_s, 1, h_Lag2Mono, N_s, h_Ustar, N_s, 0.0, h_A, N_s);
      	    // Limit the solution according to Liu
      	    Lcpu_hrl(N_s, N_E, N_F, N_G, boundaryMap, h_weight, h_monoV1D, h_J, h_A, h_Alim);
      	    // Go back to lagrange representation
      	    blasGemm("N","N", N_s, N_E*N_F, N_s, 1, h_Mono2Lag, N_s, h_Alim, N_s, 0.0, h_Ustar, N_s);
	    // Go back to conservative form
	    //Lcpu_Prim2Cons(N_s, N_E, N_F, h_Ustar, multifluid, passive, model, gamma0);
      	  }
      	} // end limiting
      }

      
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
	for(int kk = 0; kk < 2*N_F*M_T; kk++){ h_UintegF[kk] = h_UF[kk];}
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
    
    if(cpu){
      if (limiter==1){
	// Go from conservative to primitive space
	//Lcpu_Cons2Prim(N_s, N_E, N_F, h_U, multifluid, passive, model, gamma0);
	
	// // Go from nodal to modal representation
	// blasGemm("N","N", N_s, N_E*N_F, N_s, 1, h_Nod2Mod, N_s, h_U, N_s, 0.0, h_UMod, N_s);

	// // Limit the solution according to Krivodonova
	// if (blas==1) {blasCopy(N_F*N_s*N_E, h_UMod, 1, h_UModNew, 1);}    
	// else Lcpu_equal(N_s, N_E, N_F, h_UModNew, h_UMod);     // make UModNew = UMod;
	// Lcpu_hsl(N_s, N_E, N_F, boundaryMap, h_UMod, h_UModNew);
	
	// // Go back to nodal representation (slightly distors the solution at 1e-12)
	// blasGemm("N","N", N_s, N_E*N_F, N_s, 1, h_Mod2Nod, N_s, h_UModNew, N_s, 0.0, h_U, N_s);

	// Go back to conservative form
	//Lcpu_Prim2Cons(N_s, N_E, N_F, h_U, multifluid, passive, model, gamma0);
      }
      if (limiter==2){ //HR limiting
	// Go from conservative to primitive space
	//Lcpu_Cons2Prim(N_s, N_E, N_F, h_U, multifluid, passive, model, gamma0);
	// Go from lagrange to monomial representation
	blasGemm("N","N", N_s, N_E*N_F, N_s, 1, h_Lag2Mono, N_s, h_U, N_s, 0.0, h_A, N_s);
	// Limit the solution according to Liu
	Lcpu_hrl(N_s, N_E, N_F, N_G, boundaryMap, h_weight, h_monoV1D, h_J, h_A, h_Alim);
	// Go back to lagrange representation
	blasGemm("N","N", N_s, N_E*N_F, N_s, 1, h_Mono2Lag, N_s, h_Alim, N_s, 0.0, h_U, N_s);
	// Go back to conservative form
	//Lcpu_Prim2Cons(N_s, N_E, N_F, h_U, multifluid, passive, model, gamma0);
      }
    } // end limiting
    

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
  // Error calcuations
  //
  //////////////////////////////////////////////////////////////////////////

  // Initial condition
  fullMatrix<scalar> Uinit(N_s, N_E*N_F);
  if(multifluid){
    if     (simplew) init_dg_simplew_multifluid(N_s, N_E, N_F, D, model, XYZNodes, Uinit);
    else if(sodtube) init_dg_sodtube_multifluid(N_s, N_E, N_F, D, model, XYZNodes, Uinit);
    else if(contact) init_dg_contact_multifluid(N_s, N_E, N_F, D, model, XYZNodes, Uinit);
    else if(matfrnt) init_dg_matfrnt_multifluid(N_s, N_E, N_F, D, model, XYZNodes, Uinit);
    else if(sinegam) init_dg_sinegam_multifluid(N_s, N_E, N_F, D, model, XYZNodes, Uinit);
    else if(expogam) init_dg_expogam_multifluid(N_s, N_E, N_F, D, model, XYZNodes, Uinit);
  }
  else if (passive){
    if (sinephi) init_dg_sinephi_passive(N_s, N_E, N_F, D, gamma0, XYZNodes, Uinit);
  }
  scalar* h_Uinit = new scalar[N_s*N_E*N_F];  makeZero(h_Uinit,N_s*N_E*N_F);
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      for(int i = 0; i < N_s; i++){
	h_Uinit[(e*N_F+fc)*N_s+i] = Uinit(i,e*N_F+fc);}}}

  // Change to primitive variables
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      // velocity
      h_Uinit[(e*N_F+1)*N_s+i] = h_Uinit[(e*N_F+1)*N_s+i]/h_Uinit[(e*N_F+0)*N_s+i]; 
      h_U    [(e*N_F+1)*N_s+i] = h_U    [(e*N_F+1)*N_s+i]/h_U    [(e*N_F+0)*N_s+i];
      if      (multifluid){
	// gamma: get everything in terms of 1/(gamma-1)
	if      (model==0){
	  h_Uinit[(e*N_F+3)*N_s+i] = h_Uinit[(e*N_F+3)*N_s+i]/h_Uinit[(e*N_F+0)*N_s+i];
	  h_U    [(e*N_F+3)*N_s+i] = h_U    [(e*N_F+3)*N_s+i]/h_U    [(e*N_F+0)*N_s+i];}
	else if (model==1){
	  h_Uinit[(e*N_F+3)*N_s+i] = h_Uinit[(e*N_F+3)*N_s+i];
	  h_U    [(e*N_F+3)*N_s+i] = h_U    [(e*N_F+3)*N_s+i];}
	// pressure = (gamma-1)*(E-0.5 rho*v*v)
	h_Uinit[(e*N_F+2)*N_s+i] = 1.0/h_Uinit[(e*N_F+3)*N_s+i]*(h_Uinit[(e*N_F+2)*N_s+i] - 0.5*h_Uinit[(e*N_F+0)*N_s+i]*h_Uinit[(e*N_F+1)*N_s+i]*h_Uinit[(e*N_F+1)*N_s+i]); 
	h_U    [(e*N_F+2)*N_s+i] = 1.0/h_U    [(e*N_F+3)*N_s+i]*(h_U    [(e*N_F+2)*N_s+i] - 0.5*h_U    [(e*N_F+0)*N_s+i]*h_U    [(e*N_F+1)*N_s+i]*h_U    [(e*N_F+1)*N_s+i]);
      }
      else if (passive){
	// pressure = (gamma-1)*(E-0.5 rho*v*v)
	h_Uinit[(e*N_F+2)*N_s+i] = (gamma0-1)*(h_Uinit[(e*N_F+2)*N_s+i] - 0.5*h_Uinit[(e*N_F+0)*N_s+i]*h_Uinit[(e*N_F+1)*N_s+i]*h_Uinit[(e*N_F+1)*N_s+i]); 
	h_U    [(e*N_F+2)*N_s+i] = (gamma0-1)*(h_U    [(e*N_F+2)*N_s+i] - 0.5*h_U    [(e*N_F+0)*N_s+i]*h_U    [(e*N_F+1)*N_s+i]*h_U    [(e*N_F+1)*N_s+i]);
	// conservative phi
	h_Uinit[(e*N_F+3)*N_s+i] = h_Uinit[(e*N_F+3)*N_s+i]/h_Uinit[(e*N_F+0)*N_s+i]; 
	h_U    [(e*N_F+3)*N_s+i] = h_U    [(e*N_F+3)*N_s+i]/h_U    [(e*N_F+0)*N_s+i];
      }
    }
  }
    
  // Collocate the solution to the integration points
  scalar* h_Uinitg = new scalar[N_G*N_E*N_F];  makeZero(h_Uinitg,N_G*N_E*N_F);
  scalar* h_Ug     = new scalar[N_G*N_E*N_F];  makeZero(h_Ug    ,N_G*N_E*N_F);
  blasGemm("N","N", N_G, N_E*N_F, N_s, 1, h_phi, N_G, h_Uinit, N_s, 0.0, h_Uinitg, N_G);
  blasGemm("N","N", N_G, N_E*N_F, N_s, 1, h_phi, N_G, h_U    , N_s, 0.0, h_Ug    , N_G);

  // Take the cell average of the solution
  scalar* h_UinitAvg = new scalar[N_E*N_F];  makeZero(h_UinitAvg,N_E*N_F);
  scalar* h_UAvg     = new scalar[N_E*N_F];  makeZero(h_UAvg    ,N_E*N_F);
  scalar dx = XYZNodes(1,0*D+0)-XYZNodes(0,0*D+0);
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      for(int g = 0; g < N_G; g++){
	h_UinitAvg[e*N_F+fc] += h_Uinitg[(e*N_F+fc)*N_G+g]*h_J[e]*weight(g,0);
	h_UAvg    [e*N_F+fc] += h_Ug    [(e*N_F+fc)*N_G+g]*h_J[e]*weight(g,0);
      }
    }
  }

  // Go directly from nodal to modal representation
  scalar* h_UinitMod = new scalar[N_s*N_E*N_F];  makeZero(h_UinitMod,N_s*N_E*N_F);
  blasGemm("N","N", N_s, N_E*N_F, N_s, 1, h_Nod2Mod,  N_s, h_Uinit, N_s, 0.0, h_UinitMod, N_s);
  blasGemm("N","N", N_s, N_E*N_F, N_s, 1, h_Nod2Mod,  N_s, h_U    , N_s, 0.0, h_UMod    , N_s);
    
  // Get the cell average from the modal representation
  scalar* h_UinitModAvg = new scalar[N_E*N_F];  makeZero(h_UinitModAvg,N_E*N_F);
  scalar* h_UModAvg     = new scalar[N_E*N_F];  makeZero(h_UModAvg    ,N_E*N_F);
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      h_UinitModAvg[e*N_F+fc] = h_UinitMod[(e*N_F+fc)*N_s+0]*sqrt(2.0)*h_J[e]; // gamma_n = 2/(2n+1) see p45
      h_UModAvg    [e*N_F+fc] = h_UMod    [(e*N_F+fc)*N_s+0]*sqrt(2.0)*h_J[e]; // gamma_n = 2/(2n+1) see p45
    }
  }

  // Calculate the different norms of the error
  scalar E = 0;
  scalar EMod = 0;
  scalar* h_Err1      = new scalar[N_F]; makeZero(h_Err1   , N_F);
  scalar* h_Err1Mod   = new scalar[N_F]; makeZero(h_Err1Mod, N_F);
  scalar* h_Err2      = new scalar[N_F]; makeZero(h_Err2   , N_F);
  scalar* h_Err2Mod   = new scalar[N_F]; makeZero(h_Err2Mod, N_F);
  scalar* h_ErrInf    = new scalar[N_F]; makeZero(h_ErrInf   , N_F);
  scalar* h_ErrInfMod = new scalar[N_F]; makeZero(h_ErrInfMod, N_F);
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      E    = h_UinitAvg   [e*N_F+fc]-h_UAvg   [e*N_F+fc];
      EMod = h_UinitModAvg[e*N_F+fc]-h_UModAvg[e*N_F+fc];
      h_Err1[fc]    += fabs(E);
      h_Err1Mod[fc] += fabs(EMod);
      h_Err2[fc]    += E   *E;
      h_Err2Mod[fc] += EMod*EMod;
      if (h_ErrInf[fc]    < fabs(E   ))  h_ErrInf   [fc] = fabs(E);
      if (h_ErrInfMod[fc] < fabs(EMod))  h_ErrInfMod[fc] = fabs(EMod);
    }
  }
  
  // Output some stuff in a file to read by post-proc
  std::string error = "error.dat"; 
  FILE *f = fopen(error.c_str(),"w");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", h_Err1[fc]/N_E);          fprintf(f,"\n");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", h_Err1Mod[fc]/N_E);       fprintf(f,"\n");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", sqrt(h_Err2[fc]/N_E));    fprintf(f,"\n");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", sqrt(h_Err2Mod[fc]/N_E)); fprintf(f,"\n");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", h_ErrInf[fc]);            fprintf(f,"\n");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", h_ErrInfMod[fc]);         fprintf(f,"\n");
  
  // Free some stuff
  delete[] h_Uinit;
  delete[] h_Uinitg;
  delete[] h_UinitAvg;
  delete[] h_UinitModAvg;
  delete[] h_Ug;
  delete[] h_UAvg;
  delete[] h_UModAvg;
  delete[] h_UinitMod;
  delete[] h_Err1;
  delete[] h_Err1Mod;
  delete[] h_Err2;
  delete[] h_Err2Mod;
  delete[] h_ErrInf;
  delete[] h_ErrInfMod;
  


  //////////////////////////////////////////////////////////////////////////   
  //
  // Free stuff on the device
  //
  //////////////////////////////////////////////////////////////////////////   
#ifdef USE_GPU
  status = cublasShutdown();
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
#endif
  

  
  //////////////////////////////////////////////////////////////////////////   
  //
  // Free stuff on the host
  //
  //////////////////////////////////////////////////////////////////////////   

  delete[] h_Minv;
  delete[] h_phi;
  delete[] h_phiGL;
  delete[] h_phiGLinv;
  delete[] h_phi_w;
  delete[] h_dphi;
  delete[] h_dphi_w;
  delete[] h_V1D;
  delete[] h_V1Dinv;  
  delete[] h_monoV1D;
  delete[] h_monoV1Dinv;  
  delete[] h_J;
  delete[] h_invJac;
  delete[] h_U;
  delete[] h_UMod;
  delete[] h_UModNew;
  delete[] h_A;
  delete[] h_Alim;
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

  delete[] h_Nod2Mod;
  delete[] h_Mod2Nod;

  delete[] h_Lag2Mono;
  delete[] h_Mono2Lag;

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


void vandermonde1d(const int order, const fullMatrix<scalar> r, fullMatrix<scalar> &V1D){

  // Purpose : Initialize the 1D Vandermonde Matrix, V_{ij} = phi_j(r_i);
  
  V1D.resize(r.size1(),order+1);
  fullMatrix<scalar> P;
  for(int j=0;j<order+1;j++){
    JacobiP(r, 0, 0, j, P);
    for(int i=0;i<P.size1();i++) V1D(i,j) = P(i,0);
  }
}

void monovandermonde1d(const int order, const fullMatrix<scalar> r, fullMatrix<scalar> &V1D){

  // Purpose : Initialize the 1D Vandermonde Matrix, V_{ij} = (r_i)^j/factorial(j);
  
  V1D.resize(r.size1(),order+1);
  for(int j=0;j<order+1;j++){
    for(int i=0;i<r.size1();i++){
      V1D(i,j) = pow(r(i,0),j)/(scalar)factorial(j);
    }
  }
}
