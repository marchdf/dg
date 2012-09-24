#include <string>
#include <iomanip>
#include <stdio.h>
#include "stdlib.h"
#ifdef USE_GPU
#include <cutil.h>
#include <cutil_inline.h>
#include <cublas.h>
#endif
#include <time.h>
#include "fullMatrix.h"
#include "polynomialBasis.h"
#include "polynomialsJacobi.h"
#include "quadratures/Gauss.h"
#include "GmshDefines.h"
#include "simpleMesh.h"
#include <scalar_def.h>
#include <dg_functions.h>
#include <deck.h>
#include <init_cond.h>
#include <rk.h>
#include <misc.h>
#include <limiting.h>
#include <dg_solver.h>

//
// Function prototypes
//
void copyMatrixToPointer(fullMatrix<scalar> &A, scalar* h_A);
void copyMatrixToPointer(const fullMatrix<scalar> &A, scalar* h_A);
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
void average_cell_p0(const int N_s, const int N_E, const int N_F, fullMatrix<scalar> &U);

void vandermonde1d(const int order, const fullMatrix<scalar> r, fullMatrix<scalar> &V1D);
void monovandermonde1d(const int order, const fullMatrix<scalar> r, fullMatrix<scalar> &V1D);

int factorial(int n){ return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}

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
  int limiterMethod = 0;
  if      (inputs.getLimiter() == "hrl")   {limiterMethod = 1; printf("Using HR limiting\n");}
  else if (inputs.getLimiter() == "hsl")   {limiterMethod = 2; printf("Using HS limiting\n");}
  else{limiterMethod = 0; printf("No limiting\n");}

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

  // Go from lagrange to monomial basis 
  fullMatrix<scalar> Lag2Mono(N_s,N_s);
  fullMatrix<scalar> Mono2Lag(N_s,N_s);
  Lag2Mono.gemm(monoV1Dinv, phi);   // Calculate the complete nodal to modal transform = V1Dinv*phiGL
  Lag2Mono.invert(Mono2Lag);
  //Mono2Lag.gemm(phiinv,monoV1D);
  

  // //////////////////////////////////////////////////////////////////////////   
  // //
  // // Modal basis functions and change of basis
  // //
  // //////////////////////////////////////////////////////////////////////////

  // // Get the Legendre-Gauss-Lobatto points
  // fullMatrix<scalar> pointsGL;
  // JacobiGL(0,0,order,pointsGL);
  // int N_GL = pointsGL.size1();
  
  // // Evaluate the Lagrange polynomials at the GL points (needed for collocation)
  // fullMatrix<scalar> phiGL   (N_GL,N_s);
  // fullMatrix<scalar> phiGLinv(N_s,N_GL); 
  // fullMatrix<double> phiDGL  (N_GL,N_s); 
  // basis->f (pointsGL, phiDGL);
  // for(int gl = 0; gl < N_GL; gl++){
  //   for(int i = 0; i < N_s; i++){
  //     phiGL(gl,i) = (scalar)phiDGL(gl,i);
  //   }	  
  // }   
  // phiGL.invert(phiGLinv);
  
  // // Vandermonde matrix to go from modal to nodal representation
  // // u_nodal = V1D    u_modal
  // // u_modal = V1Dinv u_nodal
  // fullMatrix<scalar> V1D;
  // fullMatrix<scalar> V1Dinv;
  // vandermonde1d(order, pointsGL, V1D);
  // V1D.invert(V1Dinv);
  
  // // Go from nodal to modal
  // fullMatrix<scalar> Nod2Mod(N_s,N_s);
  // fullMatrix<scalar> Mod2Nod(N_s,N_s);
  // Nod2Mod.gemm(V1Dinv, phiGL);   // Calculate the complete nodal to modal transform = V1Dinv*phiGL
  // Nod2Mod.invert(Mod2Nod);
  
  
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
  // Setup the limiter
  //
  //////////////////////////////////////////////////////////////////////////
  scalar* h_weight  = new scalar[N_G]; makeZero(h_weight,N_G); for(int g=0; g<N_G; g++) h_weight[g] = (scalar)weight(g,0);  
  Limiting Limiter = Limiting(limiterMethod, N_s, N_E, N_F, N_G, boundaryMap, Lag2Mono, Mono2Lag, monoV1D, h_weight);
  
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
  // Initialize some stuff on the host
  //
  //////////////////////////////////////////////////////////////////////////   

  //
  //  We need to transform the data in fullMatrix to a pointer form to
  //  transfer to GPU note: it's got to be column major sorted
  //
  scalar* h_phi     = new scalar[N_G*N_s];          makeZero(h_phi,N_G*N_s);
  scalar* h_phi_w   = new scalar[N_G*N_s];          makeZero(h_phi_w,N_G*N_s);          
  scalar* h_dphi    = new scalar[D*N_G*N_s];	    makeZero(h_dphi,D*N_G*N_s);	 
  scalar* h_dphi_w  = new scalar[D*N_G*N_s];	    makeZero(h_dphi_w,D*N_G*N_s);
  scalar* h_J       = new scalar[N_E];              makeZero(h_J,N_E);                                 // not same as J!!
  scalar* h_invJac  = new scalar[N_G*D*N_E*D];      makeZero(h_invJac,N_G*D*N_E*D);                    // not same as invJac!!
  scalar* h_U       = new scalar[N_s*N_E*N_F];	    makeZero(h_U,N_s*N_E*N_F);

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
 
  printf("==== Now RK 4 steps =====\n");
  DG_SOLVER dgsolver = DG_SOLVER(D, N_F, N_E, N_s, N_G, M_T, M_s, M_G,
  				 h_phi, h_dphi, h_phi_w, h_dphi_w, h_J, h_invJac,
  				 boundaryMap, flux, model, gamma0, blas, multifluid, passive);
  RK rk4 = RK(4);
  rk4.RK_integration(Dt, N_t, output_factor,
  		     D, N_F, N_E, N_s, N_G, M_T, M_s,
  		     h_Minv, 
  		     h_U,
  		     Limiter, blas, order0, dgsolver, multifluid, passive, model, gamma0,
		     msh_lin, m);


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
  hostblasGemm('N','N', N_G, N_E*N_F, N_s, 1, h_phi, N_G, h_Uinit, N_s, 0.0, h_Uinitg, N_G);
  hostblasGemm('N','N', N_G, N_E*N_F, N_s, 1, h_phi, N_G, h_U    , N_s, 0.0, h_Ug    , N_G);

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

  // Calculate the different norms of the error
  scalar E = 0;
  scalar EMod = 0;
  scalar* h_Err1      = new scalar[N_F]; makeZero(h_Err1   , N_F);
  scalar* h_Err2      = new scalar[N_F]; makeZero(h_Err2   , N_F);
  scalar* h_ErrInf    = new scalar[N_F]; makeZero(h_ErrInf   , N_F);
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      E    = h_UinitAvg   [e*N_F+fc]-h_UAvg   [e*N_F+fc];
      h_Err1[fc]    += fabs(E);
      h_Err2[fc]    += E   *E;
      if (h_ErrInf[fc]    < fabs(E   ))  h_ErrInf   [fc] = fabs(E);
    }
  }
  
  // Output some stuff in a file to read by post-proc
  std::string error = "error.dat"; 
  FILE *f = fopen(error.c_str(),"w");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", h_Err1[fc]/N_E);          fprintf(f,"\n");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", sqrt(h_Err2[fc]/N_E));    fprintf(f,"\n");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", h_ErrInf[fc]);            fprintf(f,"\n");
  fclose(f);
  
  // Free some stuff
  delete[] h_Uinit;
  delete[] h_Uinitg;
  delete[] h_UinitAvg;
  delete[] h_Ug;
  delete[] h_UAvg;
  delete[] h_Err1;
  delete[] h_Err2;
  delete[] h_ErrInf;


  //////////////////////////////////////////////////////////////////////////   
  //
  // Free stuff on the device
  //
  //////////////////////////////////////////////////////////////////////////   
#ifdef USE_GPU
  status = cublasShutdown();
#endif
  
  
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
  delete[] h_weight;
  delete[] h_J;
  delete[] h_invJac;
  delete[] h_U;

  
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
