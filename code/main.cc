#include <string>
#include <iomanip>
#include <stdio.h>
#include "stdlib.h"
#ifdef USE_GPU
#include <cublas.h>
#endif
#ifdef USE_MPI
#include "mpi.h"
#endif
#include <ctime>
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
void get_element_types(const int order, int &msh_qua, int &msh_tri, int &msh_lin){
  if      (order==0)  {msh_qua = MSH_QUA_4;    msh_tri = MSH_TRI_3;    msh_lin = MSH_LIN_2;  }
  else if (order==1)  {msh_qua = MSH_QUA_4;    msh_tri = MSH_TRI_3;    msh_lin = MSH_LIN_2;  }
  else if (order==2)  {msh_qua = MSH_QUA_9;    msh_tri = MSH_TRI_6;    msh_lin = MSH_LIN_3;  }
  else if (order==3)  {msh_qua = MSH_QUA_16;   msh_tri = MSH_TRI_10;   msh_lin = MSH_LIN_4;  }
  else if (order==4)  {msh_qua = MSH_QUA_25;   msh_tri = MSH_TRI_15;   msh_lin = MSH_LIN_5;  }
  else if (order==5)  {msh_qua = MSH_QUA_36;   msh_tri = MSH_TRI_21;   msh_lin = MSH_LIN_6;  }
  else if (order==6)  {msh_qua = MSH_QUA_49;   msh_tri = MSH_TRI_28;   msh_lin = MSH_LIN_7;  }
  else if (order==7)  {msh_qua = MSH_QUA_64;   msh_tri = MSH_TRI_36;   msh_lin = MSH_LIN_8;  }
  else if (order==8)  {msh_qua = MSH_QUA_81;   msh_tri = MSH_TRI_45;   msh_lin = MSH_LIN_9;  }
  else if (order==9)  {msh_qua = MSH_QUA_100;  msh_tri = MSH_TRI_55;   msh_lin = MSH_LIN_10; }
  else if (order==10) {msh_qua = MSH_QUA_121;  msh_tri = MSH_TRI_66;   msh_lin = MSH_LIN_11; }
  else {printf("Invalid order number.");}
}
void average_cell_p0(const int N_s, const int N_E, fullMatrix<scalar> &U);


void vandermonde1d(const int order, const fullMatrix<scalar> r, fullMatrix<scalar> &V1D);
void monovandermonde1d(const int order, const fullMatrix<double> r, fullMatrix<scalar> &V1D);
void monovandermonde2d(const int order, const fullMatrix<double> r, fullMatrix<scalar> &V2D);
void LagMono2DTransforms(const int N_E, const int N_s, const int order, const int L2Msize1, const int L2Msize2, std::string ElemType, const fullMatrix<scalar> XYZNodes, const fullMatrix<scalar> XYZCen, fullMatrix<scalar> &Lag2Mono, fullMatrix<scalar> &Mono2Lag);
void getPowersXYZG(const int N_E, const int N_s, const int N_G, const int N_N, const int M_B, const int order, const fullMatrix<scalar> XYZG, const fullMatrix<scalar> XYZCen, const int* neighbors, const fullMatrix<scalar> shifts, scalar* powers);

int getTaylorDerIdx2DLength(const int order);
void getTaylorDerIdx2D(const int order, int* TaylorDxIdx, int* TaylorDyIdx);

void cartesian_permutations(const int order, const fullMatrix<scalar> XYZNodes, fullMatrix<scalar> &Px, fullMatrix<scalar> &Py);
void LagMono2DTransformsCartesian(const int order, const int msh_lin, const fullMatrix<scalar> Px, const fullMatrix<scalar> Py, fullMatrix<scalar> &Lag2MonoX, fullMatrix<scalar> &MonoX2MonoY, fullMatrix<scalar> &MonoY2Lag);
  
int main (int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////////////////
  //
  // Timer stuff
  //
  ////////////////////////////////////////////////////////////////////////////
  std::clock_t main_start, rk_start;
  double main_time, rk_time;
  main_start = std::clock();
  
  ////////////////////////////////////////////////////////////////////////////
  //
  // Initialize MPI if you need
  //
  ////////////////////////////////////////////////////////////////////////////
  int myid = 0; int numprocs = 1;
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  if(myid==0){printf("Total number of processors=%i and I am number %i\n",numprocs,myid);}
#endif

  ////////////////////////////////////////////////////////////////////////////
  //
  // Read arguments from deck
  //
  ////////////////////////////////////////////////////////////////////////////   
  if (argc!=3)
    if(myid==0){printf("No input deck given. Defaulting to deck.inp\n");}
  std::string deckfile = "deck.inp";
  for (int i=2;i<argc;i++){
    std::string argType = argv[i-1];
    if (argType== "-d") deckfile = argv[i];
  }
  deck inputs;  
  inputs.readDeck(deckfile.c_str());

  if(myid==0){printf("%i-dimensional problem\n",D);}
  
  // Get the blas option
#ifdef HAVE_BLAS
  if(myid==0){printf("Using BLAS\n");}
#endif

  // Get the method order
  int order = inputs.getOrder();
  bool order0 = false; if (order==0) {order0 = true; order = 1;}

  // Get the flux
#ifdef RUS
  if(myid==0){printf("Using RUSANOV\n");}
#elif HLL
  if(myid==0){printf("Using HLL\n");}
#elif ROE
  if(myid==0){printf("Using ROE\n");}
#endif
  
  // Get the mesh
  std::string fileName = inputs.getMeshfile();
  
  // Setup the limiting
  int limiterMethod = 0;
  if      (inputs.getLimiter() == "hrl")   {limiterMethod = 1; if(myid==0){printf("Using HR limiting\n");}}
  else if (inputs.getLimiter() == "myl")   {limiterMethod = 2; if(myid==0){printf("Using my limiting\n");}}
  else if (inputs.getLimiter() == "m2l")   {limiterMethod = 3; if(myid==0){printf("Using m2 limiting\n");}}
  else{limiterMethod = 0; if(myid==0){printf("No limiting\n");}}

  // Setup the initial condition type
  bool simplew = false;
  bool sodtube = false;
  bool contact = false;
  bool rhotact = false;
  bool matfrnt = false;
  bool sinegam = false;
  bool expogam = false;
  bool shckint = false;
  bool multint = false;
  bool blast1d = false;
  bool simblst = false;
  bool rarecon = false;
  bool sodcirc = false;
  bool rminstb = false;
  bool rmmulti = false;
  bool khinstb = false;
  bool khblast = false;
  bool khpertu = false;
  bool blastrm = false;
  bool sinephi = false;
  bool sodmono = false;
  if      (inputs.getInitialCondition()=="simplew") simplew = true;
  else if (inputs.getInitialCondition()=="sodtube") sodtube = true;
  else if (inputs.getInitialCondition()=="contact") contact = true;
  else if (inputs.getInitialCondition()=="rhotact") rhotact = true;
  else if (inputs.getInitialCondition()=="matfrnt") matfrnt = true;
  else if (inputs.getInitialCondition()=="sinegam") sinegam = true;
  else if (inputs.getInitialCondition()=="expogam") expogam = true;
  else if (inputs.getInitialCondition()=="shckint") shckint = true;
  else if (inputs.getInitialCondition()=="multint") multint = true;
  else if (inputs.getInitialCondition()=="blast1d") blast1d = true;
  else if (inputs.getInitialCondition()=="simblst") simblst = true;
  else if (inputs.getInitialCondition()=="rarecon") rarecon = true;
  else if (inputs.getInitialCondition()=="sodcirc") sodcirc = true;
  else if (inputs.getInitialCondition()=="rminstb") rminstb = true;
  else if (inputs.getInitialCondition()=="rmmulti") rmmulti = true;
  else if (inputs.getInitialCondition()=="khinstb") khinstb = true;
  else if (inputs.getInitialCondition()=="khblast") khblast = true;
  else if (inputs.getInitialCondition()=="khpertu") khpertu = true;
  else if (inputs.getInitialCondition()=="blastrm") blastrm = true;
  else if (inputs.getInitialCondition()=="sinephi") sinephi = true;
  else if (inputs.getInitialCondition()=="sodmono") sodmono = true;
  else{printf("Invalid initial condition setup. Correct the deck.\n");}

  // setup the boundary condition type
  int boundaryType = 0;
  if      (inputs.getBoundaryCondition()=="periodic") boundaryType = 0;
  else if (inputs.getBoundaryCondition()=="farfield") boundaryType = 1;
  else if (inputs.getBoundaryCondition()=="rflctive") boundaryType = 2;
  else if (inputs.getBoundaryCondition()=="whatever") boundaryType = 1000; // 2D doesn't need a boundary type
  else{ printf("Invalid boundary condition setup. Correct the deck.\n");}    

  //==========================================================================
  //
  //   CPU calculations
  //
  //==========================================================================

  ////////////////////////////////////////////////////////////////////////////
  //
  // Load the mesh, node coordinates, elements, interfaces, normals
  //
  ////////////////////////////////////////////////////////////////////////////   
  simpleMesh m(myid,numprocs);
  m.load(fileName.c_str());
  int elem_type;
  int face_type;
  int msh_qua;
  int msh_tri;
  int msh_lin;
  int nsides; // this offsets j in buildInterfaces function
  int N_N;    // number of neighbors to an element
  scalar refArea; // area of reference element
  get_element_types(order, msh_qua, msh_tri, msh_lin);
  if     (inputs.getElemType() == "lin"){face_type = MSH_PNT, elem_type = msh_lin; nsides = 0; N_N = 2;}
  else if(inputs.getElemType() == "tri"){face_type = msh_lin, elem_type = msh_tri; nsides = 3; N_N = 3; refArea = 0.5;}
  else if(inputs.getElemType() == "qua"){face_type = msh_lin, elem_type = msh_qua; nsides = 4; N_N = 4; refArea = 4;} 
  else printf("Invalid element type in deck");

  // Get the nodes, elements, interfaces, normals
  const fullMatrix<double> &nodes = m.getNodes();
  const std::vector<simpleElement> &elements = m.getElements(elem_type);
  m.buildInterfaces(face_type, elem_type,nsides);
  const std::vector<simpleInterface> &interfaces = m.getInterfaces();

  m.buildNormals(face_type, elem_type);
  const fullMatrix<scalar> &normals = m.getNormals();

  m.buildElementMap(elem_type);
  const std::map<int,int> &ElementMap = m.getElementMap();

  m.buildCommunicators(elem_type); // build the indexes to map the ghost elements to my partition
  const std::map<int,int> & ghostElementMap = m.getGhostElementMap();
  int* h_ghostInterfaces  = m.getGhostInterfaces();
  int* h_ghostElementSend = m.getGhostElementSend();
  int* h_ghostElementRecv = m.getGhostElementRecv();

 
  ////////////////////////////////////////////////////////////////////////////   
  //
  // Generer les fonctions de formes, get integration points, weights
  //
  ////////////////////////////////////////////////////////////////////////////   

  // Elements
  const polynomialBasis *basis  = polynomialBases::find(elem_type);  // for the element
  const std::vector<std::vector<int> > &closures = basis->closures;
  fullMatrix<double> points, weight;
  if     (inputs.getElemType() == "lin") gaussIntegration::getLine(order*2+1, points, weight);
  else if(inputs.getElemType() == "tri") gaussIntegration::getTriangle(order*2+1, points, weight);
  else if(inputs.getElemType() == "qua") gaussIntegration::getQuad(order*2+1, points, weight);

  // Faces
  const polynomialBasis *basisF; // for the edges
#ifdef TWOD
  basisF = polynomialBases::find (msh_lin);
#endif
  fullMatrix<double> pointsF, weightF;
#ifdef ONED
  pointsF.resize(1,3); weightF.resize(1,1); weightF(0,0)=1;
#elif TWOD
  gaussIntegration::getLine(order*2+1, pointsF, weightF);
#endif

  //////////////////////////////////////////////////////////////////////////   
  //
  // Define some numbers for clarity
  //
  //////////////////////////////////////////////////////////////////////////   

  int N_s = elements[0].getNbNodes(); // number of nodes on an element          (i index)
  int M_s = 1;                        // number of nodes on a face              (j index)
#ifdef ONED
  M_s = 1;
#elif TWOD
  M_s = order + 1;
#endif
  int N_T = basis->numFaces;          // number of faces per element            
  int N_E = elements.size();          // number of elements                     (e index)
  int M_T = interfaces.size();        // number of faces                        (t index)
  int N   = N_E*N_s;                  // number of dof for a DG
  int N_G = points.size1();           // number of integration points           (g index)
  int M_G = pointsF.size1();          // number of integration points on a face (g index)
  int N_ghosts = m.getNbGhostInterfaces();

  if(myid==0){
    if(order0) printf("order %i\n",0); else printf("order %i\n",order);
    printf("N_s %i\n",N_s);
    printf("M_s %i\n",M_s);
    printf("N_T %i\n",N_T);
    printf("N_E %i\n",N_E);
    printf("M_T %i\n",M_T);
    printf("N_G %i\n",N_G);
    printf("M_G %i\n",M_G);
    printf("N_ghosts %i\n",N_ghosts);
  }
  
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
  
  // Faces
  fullMatrix<scalar> psi(M_G, M_s);
  fullMatrix<scalar> dpsi(M_G*DF,M_s);
#ifdef ONED
  psi(0,0)  = 1; dpsi(0,0) = 0;
#elif TWOD
  fullMatrix<double> psiD(M_G, M_s);
  basisF->f (pointsF,psiD);
  for(int g = 0; g < M_G; g++){
    for(int j = 0; j < M_s; j++){
      psi(g,j) = (scalar)psiD(g,j);
    }
  }
  double gradsF[M_s][3];
  for(int g = 0; g < M_G; g++){
    basisF->df(pointsF(g,0),pointsF(g,1),pointsF(g,2),gradsF);
    for(int alpha = 0; alpha < DF; alpha ++){
      for(int j = 0; j < M_s; j++){
	dpsi(g*DF+alpha,j) = (scalar)gradsF[j][alpha];  // see paper for indexing p.6
      }	  
    }    
  }
#endif


  
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

  // Faces
  fullMatrix<scalar> psi_w (M_G,M_s);
  for(int g = 0; g < M_G; g++){
    for(int j = 0; j < M_s; j++){
      psi_w(g,j) = (scalar)psi(g,j)*weightF(g,0);
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

  // Element centroids
  fullMatrix<scalar> XYZCen = m.getElementCentroids(N_E, N_N, XYZNodes);

  // Is this a cartesian mesh?
  bool cartesian = m.iscartesian(inputs.getElemType(),elem_type);
  if(myid==0){
    if  (cartesian) printf("Structured mesh\n");
    else            printf("Unstructured mesh\n");
  }
  
  // Faces
  fullMatrix<scalar> XYZNodesF (M_s, M_T*2*D);
  fullMatrix<scalar> XYZGF (M_G, M_T*2*D);
  for (int t = 0; t < M_T; t++) {
    const simpleInterface &face = interfaces[t];
    for(int d = 0; d < 2; d++){
      const simpleElement *el = face.getElement(d);
      if(el->getPartition()==myid){
	//if(el!=NULL){
  	int id = face.getClosureId(d);
  	const std::vector<int> &cl = closures[id];
  	for(int j = 0; j < M_s; j++){
  	  for(int alpha = 0; alpha < D; alpha++){
  	    XYZNodesF(j,(t*2+d)*D+alpha) = XYZNodes(cl[j],ElementMap.at(el->getId())*D+alpha);
  	  }
  	}
      }
      else{
	for(int j = 0; j < M_s; j++){
  	  for(int alpha = 0; alpha < D; alpha++){
	    XYZNodesF(j,(t*2+1)*D+alpha) = XYZNodesF(j,(t*2+0)*D+alpha);
	  }
	}
      }
    }
  }
  XYZGF.gemm (psi, XYZNodesF);

  
  //////////////////////////////////////////////////////////////////////////   
  //
  // Build the boundary map (must be done after the normals)
  //
  //////////////////////////////////////////////////////////////////////////
  m.buildBoundary();
  int* h_boundaryMap = m.getBoundaryMap();
  int* h_boundaryIdx = m.getBoundaryIdx();
  int M_B = m.getBoundarySize();
  if(myid==0){printf("M_B %i\n",M_B);}
  
  //////////////////////////////////////////////////////////////////////////   
  //
  // Build neighbors map (must be done after ElementMap and boundaryMap)
  //
  //////////////////////////////////////////////////////////////////////////

  // Get the coordinate shifts to bring the neighbors of the elements
  // on the boundaries to the correct (x,y) location. This is used
  // when evaluating polynomials in neighboring elements.
// #ifdef ONED
//   m.buildBoundaryElementShift1D(N_s, N_E, XYZNodes);
// #elif TWOD
//   m.buildBoundaryElementShift2D(order, XYZNodesF, ElementMap);
// #endif
//   fullMatrix<scalar> shifts = m.getShifts();

  m.buildNeighbors(N_N, N_E);
  int* h_neighbors = m.getNeighbors();

  //////////////////////////////////////////////////////////////////////////   
  //
  // Monomial to Lagrange basis transforms
  //
  //////////////////////////////////////////////////////////////////////////
#ifdef ONED
  fullMatrix<scalar> monoV;
  fullMatrix<scalar> monoVinv;
  monovandermonde1d(order, points, monoV);
  monoV.invert(monoVinv);

  // Go from lagrange to monomial basis (in ref space)
  fullMatrix<scalar> Lag2Mono(N_s,N_s);
  fullMatrix<scalar> Mono2Lag(N_s,N_s);
  Lag2Mono.gemm(monoVinv, phi);   // Calculate the complete nodal to modal transform = V1Dinv*phiGL
  Lag2Mono.invert(Mono2Lag);

#elif TWOD
  //
  // Structured/uniform mesh
  //
  fullMatrix<scalar> Lag2MonoX;
  fullMatrix<scalar> MonoX2MonoY;
  fullMatrix<scalar> MonoY2Lag;
  fullMatrix<scalar> monoV;
  if(cartesian){
    fullMatrix<scalar> Px;
    fullMatrix<scalar> Py;
    cartesian_permutations(order,XYZNodes, Px,Py);
    LagMono2DTransformsCartesian(order, msh_lin, Px, Py, Lag2MonoX, MonoX2MonoY, MonoY2Lag);
    monovandermonde1d(order, pointsF, monoV);
  }

  // //
  // // Unstructured mesh limiting
  // //
  // int L2Msize1, L2Msize2; // number of rows and columns in Lag2Mono transforms
  // if     (inputs.getElemType()== "tri"){L2Msize1 = N_s; L2Msize2 = N_s;}
  // else if(inputs.getElemType()== "qua"){L2Msize1 = (2*order+1)*(order+1); L2Msize2 = N_s;}
  // fullMatrix<scalar> Lag2Mono(N_E, L2Msize1*L2Msize2);
  // fullMatrix<scalar> Mono2Lag(N_E, L2Msize2*L2Msize1);
  // scalar* h_powersXYZG;
  // if(!cartesian){   // Go from lagrange to monomial basis (in physical space)
  //   LagMono2DTransforms(N_E, N_s, order, L2Msize1, L2Msize2, inputs.getElemType(), XYZNodes, XYZCen, Lag2Mono, Mono2Lag);
  
  //   // Get the powers of the physical nodes and neighbors
  //   // required for limiting in 2D
  //   if     (inputs.getElemType() == "tri"){
  //     h_powersXYZG = new scalar[N_s*N_G*(N_N+1)*N_E];
  //     getPowersXYZG(N_E, N_s, N_G, N_N, M_B, order, XYZG, XYZCen, h_neighbors, shifts, h_powersXYZG);}
  //   else if(inputs.getElemType() == "qua"){
  //     h_powersXYZG = new scalar[L2Msize1*N_G*(N_N+1)*N_E];
  //     getPowersXYZG(N_E, L2Msize1, N_G, N_N, M_B, 2*order, XYZG, XYZCen, h_neighbors, shifts, h_powersXYZG);}
  // }
  
#endif
    
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
  // Initialize the unknowns
  //
  //////////////////////////////////////////////////////////////////////////
  
  fullMatrix<scalar> U(N_s, N_E*N_F);
  fullMatrix<scalar> Us(N_s, N_E*N_F);
  fullMatrix<scalar> Ustar(N_s, N_E*N_F);
#ifdef MULTIFLUID
  if     (simplew) init_dg_simplew_multifluid(N_s, N_E, XYZNodes, U);
  else if(sodtube) init_dg_sodtube_multifluid(N_s, N_E, XYZNodes, U);
  else if(sodmono) init_dg_sodmono_multifluid(N_s, N_E, XYZNodes, U);
  else if(contact) init_dg_contact_multifluid(N_s, N_E, XYZNodes, U);
  else if(rhotact) init_dg_rhotact_multifluid(N_s, N_E, XYZNodes, U);
  else if(matfrnt) init_dg_matfrnt_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(sinegam) init_dg_sinegam_multifluid(N_s, N_E, XYZNodes, U);
  else if(expogam) init_dg_expogam_multifluid(N_s, N_E, XYZNodes, U);
  else if(shckint) init_dg_shckint_multifluid(N_s, N_E, XYZNodes, U);
  else if(multint) init_dg_multint_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(blast1d) init_dg_blast1d_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(simblst) init_dg_simblst_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(rarecon) init_dg_rarecon_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(sodcirc) init_dg_sodcirc_multifluid(N_s, N_E, XYZNodes, U);
  else if(rminstb) init_dg_rminstb_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(rmmulti) init_dg_rmmulti_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(khinstb) init_dg_khinstb_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(khblast) init_dg_khblast_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(khpertu) init_dg_khpertu_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(blastrm) init_dg_blastrm_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
#elif PASSIVE
  if (sinephi) init_dg_sinephi_passive(N_s, N_E, XYZNodes, U);
  if (sodmono) init_dg_sodmono_passive(N_s, N_E, XYZNodes, U);
#endif

  if (order0) average_cell_p0(N_s, N_E, U);
  

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
  dg_jacobians_elements(N_G, N_E, XYZNodes, dphi, Jac, invJac, J, invJ);
  
  // Faces
  fullMatrix<scalar> JacF(M_T*2*D,M_G*DF);
  fullMatrix<scalar> JF(M_T*2,1);            // determinant of the Jacobian
  fullMatrix<scalar> invJF(M_T*2,1);         // determinant of the inverse Jacobian
  dg_jacobians_face(M_T, XYZNodesF, dpsi, JacF, JF, invJF);
  
  //////////////////////////////////////////////////////////////////////////   
  // 
  // Calculate the inverse mass matrices
  //
  //////////////////////////////////////////////////////////////////////////

  scalar* h_Minv = new scalar[N_s*N_s*N_E];
  dg_inverse_mass_matrix(order, elem_type, inputs.getElemType(), N_s, N_E, XYZNodes, h_Minv);

  //////////////////////////////////////////////////////////////////////////   
  // 
  // Build the map
  //
  //////////////////////////////////////////////////////////////////////////
  int* h_map = new int[M_s*M_T*N_F*2];
  int* h_invmap = new int[M_s*N_N*N_E*N_F*2];
  dg_mappings(myid, M_s, M_T, N_s, N_E, N_N, interfaces, ElementMap, closures, h_map, h_invmap);
  
  //==========================================================================
  //
  //   GPU calculations
  //
  //==========================================================================
#ifdef USE_GPU
  // Choose the device
  //cudaSetDevice(0);

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

  scalar* h_psi     = new scalar[M_G*M_s];	    makeZero(h_psi,M_G*M_s);	 
  scalar* h_psi_w   = new scalar[M_G*M_s];	    makeZero(h_psi_w,M_G*M_s);	 
  scalar* h_J       = new scalar[N_E];              makeZero(h_J,N_E);               // not same as J!!
  scalar* h_invJac  = new scalar[N_G*D*N_E*D];      makeZero(h_invJac,N_G*D*N_E*D);  // not same as invJac!!
  scalar* h_JF      = new scalar[2*M_T];            makeZero(h_JF, D*M_T);
  scalar* h_normals = new scalar[D*M_T];	    makeZero(h_normals,D*M_T);	 
  scalar* h_U       = new scalar[N_s*N_E*N_F];	    makeZero(h_U,N_s*N_E*N_F);

  // copy from the fullMatrix to the pointer format (column major)
  phi.copyMatrixToPointer(h_phi);
  phi_w.copyMatrixToPointer(h_phi_w);
  dphi.copyMatrixToPointer(h_dphi);
  dphi_w.copyMatrixToPointer(h_dphi_w);
  psi.copyMatrixToPointer(h_psi);
  psi_w.copyMatrixToPointer(h_psi_w);
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
  for(int t = 0; t < M_T; t++){
    for(int d = 0; d < 2; d++){
      h_JF[t*2+d] = JF(t*2+d,0);
    }
  }
  normals.copyMatrixToPointer(h_normals);
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      for(int i = 0; i < N_s; i++){
  	h_U[(e*N_F+fc)*N_s+i]= U(i,e*N_F+fc);
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////   
  //
  // Setup the limiter
  //
  //////////////////////////////////////////////////////////////////////////
  scalar* h_weight  = new scalar[N_G]; makeZero(h_weight,N_G); for(int g=0; g<N_G; g++) h_weight[g] = (scalar)weight(g,0);  
#ifdef ONED
  Limiting Limiter = Limiting(limiterMethod, N_s, N_E, N_G, N_N, h_neighbors, Lag2Mono, Mono2Lag, monoV, h_weight);
#elif TWOD

  //
  // Structured mesh
  //
  //if(cartesian){
  scalar* h_weightF  = new scalar[M_G]; makeZero(h_weightF,M_G); for(int g=0; g<M_G; g++) h_weightF[g] = (scalar)weightF(g,0);  
  Limiting Limiter = Limiting(limiterMethod, N_s, N_E, N_G, order, cartesian, N_N, N_ghosts, h_neighbors, Lag2MonoX, MonoX2MonoY, MonoY2Lag, monoV, h_ghostElementSend, h_ghostElementRecv, h_weightF);
  delete[] h_weightF;
  //}
  
  // //
  // // Unstructured mesh
  // //
  // if(!cartesian){
  //   int L = getTaylorDerIdx2DLength(order);
  //   int* h_TaylorDxIdx = new int[L];
  //   int* h_TaylorDyIdx = new int[L];
  //   getTaylorDerIdx2D(order, h_TaylorDxIdx, h_TaylorDyIdx);
  //   Limiting Limiter = Limiting(limiterMethod, N_s, N_E, N_F, N_G, N_N, L, order, L2Msize1, L2Msize2, h_neighbors, Lag2Mono, Mono2Lag, XYZCen, h_powersXYZG, h_weight, refArea, h_TaylorDxIdx, h_TaylorDyIdx);
  //   delete[] h_TaylorDxIdx; delete[] h_TaylorDyIdx;
  //   delete[] h_powersXYZG;
  // }
 
#endif

  //////////////////////////////////////////////////////////////////////////   
  //
  // Solve the problem on the CPU/GPU.
  //
  //////////////////////////////////////////////////////////////////////////   
  double DtOut = inputs.getOutputTimeStep(); 
  double Tf    = inputs.getFinalTime();
  m.setDx(N_N,N_E,XYZCen,XYZNodes);
  scalar CFL   = inputs.getCFL()*m.getDx()/(2.0*order+1);//according to Eder. I used to have: (2.0*order+1);

  if(myid==0){printf("==== Now RK 4 steps =====\n");}
  DG_SOLVER dgsolver = DG_SOLVER(N_E, N_s, N_G, N_N, M_T, M_s, M_G, M_B, N_ghosts,
  				 h_map, h_invmap, h_ghostInterfaces, h_phi, h_dphi, h_phi_w, h_dphi_w, h_psi, h_psi_w, h_J, h_invJac, h_JF, h_weight, h_normals,
  				 h_boundaryMap, h_boundaryIdx);
  RK rk4 = RK(4);




  /* Your algorithm here */
  rk_start = std::clock();
  rk4.RK_integration(DtOut, Tf, CFL,
  		     N_E, N_s, N_G, M_T, M_s,
  		     h_Minv, 
  		     h_U,
  		     Limiter, order0, dgsolver,
  		     elem_type, m);
  rk_time = ( std::clock() - rk_start ) / (double) CLOCKS_PER_SEC;
  printf("RK time = %20.16e for proc %i\n", rk_time, myid);

  //////////////////////////////////////////////////////////////////////////   
  //
  // Error calcuations
  //
  //////////////////////////////////////////////////////////////////////////

#ifdef ERROR
  // Initial condition
  fullMatrix<scalar> Uinit(N_s, N_E*N_F);
#ifdef MULTIFLUID
  if     (simplew) init_dg_simplew_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(sodtube) init_dg_sodtube_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(contact) init_dg_contact_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(rhotact) init_dg_rhotact_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(matfrnt) init_dg_matfrnt_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(sinegam) init_dg_sinegam_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(expogam) init_dg_expogam_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(shckint) init_dg_shckint_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(rarecon) init_dg_rarecon_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(sodcirc) init_dg_sodcirc_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(rminstb) init_dg_rminstb_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(rmmulti) init_dg_rmmulti_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(khblast) init_dg_khblast_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(khpertu) init_dg_khpertu_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(blastrm) init_dg_blastrm_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
#elif PASSIVE
  if (sinephi) init_dg_sinephi_passive(N_s, N_E, XYZNodes, Uinit);
  if (sodmono) init_dg_sodmono_passive(N_s, N_E, XYZNodes, Uinit);
#endif
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
#ifdef MULTIFLUID
      // gamma: get everything in terms of 1/(gamma-1)
#ifdef GAMCONS
      h_Uinit[(e*N_F+3)*N_s+i] = h_Uinit[(e*N_F+3)*N_s+i]/h_Uinit[(e*N_F+0)*N_s+i];
      h_U    [(e*N_F+3)*N_s+i] = h_U    [(e*N_F+3)*N_s+i]/h_U    [(e*N_F+0)*N_s+i];
#elif GAMNCON
      h_Uinit[(e*N_F+3)*N_s+i] = h_Uinit[(e*N_F+3)*N_s+i];
      h_U    [(e*N_F+3)*N_s+i] = h_U    [(e*N_F+3)*N_s+i];
#endif
      // pressure = (gamma-1)*(E-0.5 rho*v*v)
      h_Uinit[(e*N_F+2)*N_s+i] = 1.0/h_Uinit[(e*N_F+3)*N_s+i]*(h_Uinit[(e*N_F+2)*N_s+i] - 0.5*h_Uinit[(e*N_F+0)*N_s+i]*h_Uinit[(e*N_F+1)*N_s+i]*h_Uinit[(e*N_F+1)*N_s+i]); 
      h_U    [(e*N_F+2)*N_s+i] = 1.0/h_U    [(e*N_F+3)*N_s+i]*(h_U    [(e*N_F+2)*N_s+i] - 0.5*h_U    [(e*N_F+0)*N_s+i]*h_U    [(e*N_F+1)*N_s+i]*h_U    [(e*N_F+1)*N_s+i]);
#elif PASSIVE
      // pressure = (gamma-1)*(E-0.5 rho*v*v)
      scalar gamma = constants::GLOBAL_GAMMA;
      h_Uinit[(e*N_F+2)*N_s+i] = (gamma-1)*(h_Uinit[(e*N_F+2)*N_s+i] - 0.5*h_Uinit[(e*N_F+0)*N_s+i]*h_Uinit[(e*N_F+1)*N_s+i]*h_Uinit[(e*N_F+1)*N_s+i]); 
      h_U    [(e*N_F+2)*N_s+i] = (gamma-1)*(h_U    [(e*N_F+2)*N_s+i] - 0.5*h_U    [(e*N_F+0)*N_s+i]*h_U    [(e*N_F+1)*N_s+i]*h_U    [(e*N_F+1)*N_s+i]);
      // conservative phi
      h_Uinit[(e*N_F+3)*N_s+i] = h_Uinit[(e*N_F+3)*N_s+i]/h_Uinit[(e*N_F+0)*N_s+i]; 
      h_U    [(e*N_F+3)*N_s+i] = h_U    [(e*N_F+3)*N_s+i]/h_U    [(e*N_F+0)*N_s+i];
#endif
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

#endif

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
  delete[] h_ghostInterfaces;
  delete[] h_ghostElementSend;
  delete[] h_ghostElementRecv;
  delete[] h_boundaryMap;
  delete[] h_boundaryIdx;
  delete[] h_neighbors;
  delete[] h_Minv;
  delete[] h_map;
  delete[] h_invmap;
  delete[] h_phi;
  delete[] h_phi_w;
  delete[] h_dphi;
  delete[] h_dphi_w;
  delete[] h_psi;
  delete[] h_psi_w;
  delete[] h_weight;
  delete[] h_J;
  delete[] h_JF;
  delete[] h_invJac;
  delete[] h_normals;
  delete[] h_U;

  //////////////////////////////////////////////////////////////////////////   
  //
  // Finalize MPI
  //
  //////////////////////////////////////////////////////////////////////////
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
  MPI_Finalize();
#endif

  //////////////////////////////////////////////////////////////////////////   
  //
  // End timer
  //
  //////////////////////////////////////////////////////////////////////////
  main_time = ( std::clock() - main_start ) / (double) CLOCKS_PER_SEC;
  printf("Main time = %20.16e for proc %i\n", main_time, myid);

  return 0;
}// end main



//
// Function definitions
//
void average_cell_p0(const int N_s, const int N_E, fullMatrix<scalar> &U){
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

void monovandermonde1d(const int order, const fullMatrix<double> r, fullMatrix<scalar> &V1D){

  // Purpose : Initialize the 1D Vandermonde Matrix, V_{ij} = (r_i)^j/factorial(j);
  
  V1D.resize(r.size1(),order+1);
  for(int j=0;j<order+1;j++){
    for(int i=0;i<r.size1();i++){
      V1D(i,j) = pow(r(i,0),j)/(scalar)factorial(j);
    }
  }
}

void monovandermonde2d(const int order, const fullMatrix<double> r, fullMatrix<scalar> &V2D){

  // Purpose : Initialize the 2D Vandermonde Matrix, V = (x_i)^nx (y_i)^ny / (factorial(nx)*factorial(ny));
  // ith line = [1, x_i, y_i, x_i^2/2!, x_i*y_i, y_i^2/2!, x_i^3/3!, x_i^2*y_i/2!, x_i*y_i^2/2!, y_i^3/3!, ...]

  V2D.resize(r.size1(),(int)((order+1)*(order+2)/2.0));

  for(int i=0;i<r.size1();i++){
    int offset = 0;
    for(int p=0; p<order+1; p++){
      for(int k=0; k<p+1; k++){
	int nx = p-k;
	int ny =   k;
	//printf("line(%i) for p=%i: nx=%i and ny=%i (column=%i)\n",i,p,nx,ny,offset+k);
	V2D(i,offset+k) = pow(r(i,0),nx)*pow(r(i,1),ny)/(scalar)(factorial(nx)*factorial(ny));
      }
      offset = offset + p + 1;
    }
  }
}

// Returns the transforms (and inverse) from Lagrange to Taylor
// polynomials NB: In >< to the 1D transforms, these are in the physical
// space! So there is one transform per element
void LagMono2DTransforms(const int N_E, const int N_s, const int order, const int L2Msize1, const int L2Msize2, std::string ElemType, const fullMatrix<scalar> XYZNodes, const fullMatrix<scalar> XYZCen, fullMatrix<scalar> &Lag2Mono, fullMatrix<scalar> &Mono2Lag){

  fullMatrix<scalar> M2L;
  fullMatrix<scalar> L2M;
  fullMatrix<double> points(N_s,D);
  // Modifications if you are dealing with a quadrangle Calculate
  // the A matrix to go from partial taylor polynomial to full
  // taylor (T_f = A * T_p) or from full vandermonde matrix to partial
  // vandermonde matrix (V = V_f * A)
  fullMatrix<scalar> A;
  if(ElemType == "qua"){
    A.resize(L2Msize1,L2Msize2); A.scale((scalar)0.0);
    int i=0,j=0,cnt=0;
    for(int idx = 0; idx < 2*order+1; idx++){
      if(idx<=order){
	for(int k=0; k<=idx; k++){ A(i,j) = 1; i++; j++;}
      }
      else if(idx>order){
	for(int k=0; k<=idx; k++){
	  if((cnt<k)&&(k<idx-cnt)){ A(i,j) = 1; i++; j++;}
	  else{ i++;}
	}
	cnt++;
      }
    }
  }
  fullMatrix<scalar> At = A.transpose(); // transpose of A
  
  // Loop on elements
  for(int e = 0; e < N_E; e++){

    // Get the points
    for(int i = 0; i < N_s; i++)
      for(int alpha = 0; alpha < D; alpha++)
	points(i,alpha) = XYZNodes(i,e*D+alpha)-XYZCen(e,alpha);

    // Get the power matrix
    if(ElemType == "tri"){
      monovandermonde2d(order, points, M2L);
      M2L.invert(L2M);
    }
    else if(ElemType == "qua"){
      fullMatrix<scalar> V(N_s,N_s);
      fullMatrix<scalar> Vinv(N_s,N_s);
      fullMatrix<scalar> Vf; // full vandermonde matrix
      monovandermonde2d(2*order, points, Vf);
      V.gemm(Vf,A);         // strip some columns out of Vf
      V.invert(Vinv);       // get the inverse

      // build the transform matrices
      L2M.resize(L2Msize1, L2Msize2);
      M2L.resize(L2Msize2, L2Msize1);
      L2M.gemm(A,Vinv);  // T_f = (A*Vinv)*phi
      M2L.gemm(V,At);    // phi = (V*At)  *T_f
    }

    // Store them in the large transform matrix
    for(int i = 0; i < L2Msize1; i++){
      for(int j = 0; j < L2Msize2; j++){
    	Lag2Mono(e,i*L2Msize2+j) = L2M(i,j); // size: L2Msize1 x L2Msize2
	Mono2Lag(e,j*L2Msize1+i) = M2L(j,i); // size: L2Msize2 x L2Msize1
      }
    }    
  } // element loop
}


// Get the powers of XZYG-XYZCen for each element and his neighbors
// This is precalculated for increased speed in 2D limiting
void getPowersXYZG(const int N_E, const int N_s, const int N_G, const int N_N, const int M_B, const int order, const fullMatrix<scalar> XYZG, const fullMatrix<scalar> XYZCen, const int* neighbors, const fullMatrix<scalar> shifts, scalar* powers){

  fullMatrix<scalar> V;
  fullMatrix<double> points(N_G,D);
  int el = 0; // index of the neighboring element (or itself)
  scalar* XYZshift = new scalar[D]; 
  
  for(int e = 0; e < N_E; e++){
    for(int nn = 0; nn < N_N+1; nn++){
      if (nn == 0) {el = e;}
      else         {el = neighbors[e*N_N+nn-1];} // -1 bc nn starts with the current element

      // Get potentional coordinate shift for the neighbor if this
      // element is on the boundary
      bool flag = false; int bidx = 0;
      for(int b=0; b<M_B; b++)
	if ((e==shifts(b,0))&&(el==shifts(b,1))){ flag = true; bidx = b; break;}
      if   (flag) for(int alpha=0; alpha<D; alpha++) XYZshift[alpha] = shifts(bidx,2+alpha);
      else        for(int alpha=0; alpha<D; alpha++) XYZshift[alpha] = 0;
      
      // Get the points to evaluate the Taylor polynomial
      // in the element or the neighboring elements
      for(int g = 0; g < N_G; g++)
  	for(int alpha = 0; alpha < D; alpha++)
  	  points(g,alpha) = XYZG(g,el*D+alpha) - XYZCen(e,alpha) + XYZshift[alpha];

      // Get the powers of these points
      monovandermonde2d(order, points, V);

      // Store in these in the large matrix
      for(int g = 0; g < N_G; g++)
  	for(int i = 0; i < N_s; i++)
  	  powers[((e*(N_N+1)+nn)*N_G+g)*N_s+i] = V(g,i);
    } // nn loop
  } // e loop

  // Free some stuff
  delete[] XYZshift;
}


// Get the x and y derivative index. Basically gives you the index for
// various derivatives of a Taylor polynomial. order is the DG order
// dxIdx =[idx for 0th derivative wrt x, idx for 1st derivative wrt x, ...]
int getTaylorDerIdx2DLength(const int order){
  // Calculate the length of these indexes
  int L = 0;
  for(int p = order; p >= 0; p--) L+=(p+1)*(p+2)/2;  // DOES THIS WORK ONLY FOR TRIANGLES?
  return L;
}

void getTaylorDerIdx2D(const int order, int* TaylorDxIdx, int* TaylorDyIdx){

  // Fill these vectors with the appropriate indexes
  // wrt x
  int countx = 0;
  int kstart = 0;
  for(int pstart = 0; pstart <= order; pstart++){ // loop on all derivative orders
    int offset = pstart*(pstart+1)/2;
    int kend = 0;
    for(int p = pstart; p <= order; p++){ 
      for(int k = kstart; k <= kend+kstart; k++){
  	TaylorDxIdx[countx] = offset + k;
  	countx++;
      }
      offset = offset + p + 1;
      kend++;
    }
  }// end pstart loop

  // wrt y
  int county = 0;
  kstart = 0;
  for(int pstart = 0; pstart <= order; pstart++){ // loop on all derivative orders
    int offset = pstart*(pstart+1)/2;
    int kend = 0;
    for(int p = pstart; p <= order; p++){ 
      for(int k = kstart; k <= kend+kstart; k++){
  	TaylorDyIdx[county] = offset + k;
  	county++;
      }
      offset = offset + p + 1;
      kend++; 
    }
    kstart++;
  }// end pstart loop
}


void cartesian_permutations(const int order, const fullMatrix<scalar> XYZNodes, fullMatrix<scalar> &Px, fullMatrix<scalar> &Py){

  /* This function should be used for a cartesian mesh.  The idea is
   to find the permutation matrix to go from the numbering system in
   gmsh to a system where the nodes are numbered in increasing order
   for increasing x and decreasing y. Best shown by example: for p=2,
   the values of U are stored at the following indexes:

    0--4--1                                     3--5--4
    |  |  |                                     |  |  |
    7--8--5   but we want (for 2D limiting) =>  6--8--7
    |  |  |                                     |  |  |
    3--6--2                                     0--2--1

    Therefore:
        U(ordered in x) = Px*U(original)
	[3 2 6 0 1 4 7 5 8]' = Px * [0 1 2 3 4 5 6 7 8]

	U(ordered in y) = Py*U(original)
	[3 0 7 2 1 5 6 4 8]' = Px * [0 1 2 3 4 5 6 7 8]

    And the inverse transform is given by Px' and Py'
    (bc inverse of permutation matrix is the transpose)
  */

  // Allocate some resources
  int nlvl = order+1;
  int N_s = nlvl*nlvl;
  Px.resize(N_s,N_s);
  Py.resize(N_s,N_s);
  int* tmp = new int[nlvl];
  int pcnt = 0;

  fullMatrix<scalar> xy(N_s,2);
  for(int i=0;i<N_s;i++){
    xy(i,0) = XYZNodes(i,0*2+0); // use the first element (for example...)
    xy(i,1) = XYZNodes(i,0*2+1);
  }

  //
  // Initialize the levels of constant x and y
  //
  double* xlvl = new double[nlvl];
  double* ylvl = new double[nlvl]; 

  // Find xmax, xmin, ymax, ymin
  scalar xmax=xy(0,0),xmin=xy(0,0),ymin=xy(0,1),ymax=xy(0,1);
  for(int i=1; i<N_s; i++){
    xmax = MAX(xmax,xy(i,0));
    xmin = MIN(xmin,xy(i,0));
    ymax = MAX(ymax,xy(i,1));
    ymin = MIN(ymin,xy(i,1));
  }

  // Determine the spacing of the lvls of constant x and y
  scalar dx = (xmax-xmin)/order;
  scalar dy = (ymax-ymin)/order;
  
  // Build the levels of constant x and y. First do the min, then the
  // max, then fill in the rest in increasing order.
  xlvl[0] = xmin; xlvl[1] = xmax;
  ylvl[0] = ymin; ylvl[1] = ymax;
  for(int p=2;p<nlvl;p++){
    xlvl[p] = xmin+(p-1)*dx;
    ylvl[p] = ymin+(p-1)*dy;
  }

  //
  // Get the permutations wrt x
  // 
  pcnt=0;
  // Loop on all the y levels
  for(int ycnt=0; ycnt<nlvl; ycnt++){
    int cnt = 0;
    // Loop to find all the nodes at that ylvl
    for(int i=0; i<N_s; i++){
      if(fabs(ylvl[ycnt]-xy(i,1))<1e-9) {tmp[cnt] = i; cnt++;}
    }// end loop on nodes

    // Sort the nodes at this level in ascending x coord
    for(int xcnt=0; xcnt<nlvl; xcnt++){
      for(int i=0;i<nlvl;i++){
	if(fabs(xlvl[xcnt]-xy(tmp[i],0)) < 1e-9){
	  Px(pcnt,tmp[i]) = 1;
	  pcnt++;
	}
      }
    }    
  } // end ycnt loop

  //
  // Get permutations wrt y
  //
  pcnt = 0;
  // Loop on all the x levels
  for(int xcnt=0; xcnt<nlvl; xcnt++){
    int cnt = 0;
    // Loop to find all the nodes at that xlvl
    for(int i=0; i<N_s; i++){
      if(fabs(xlvl[xcnt]-xy(i,0))<1e-9) {tmp[cnt] = i; cnt++;}
    }// end loop on nodes

    // Sort the nodes at this level in ascending y coord
    for(int ycnt=0; ycnt<nlvl; ycnt++){
      for(int i=0;i<order+1;i++){
	if(fabs(ylvl[ycnt]-xy(tmp[i],1)) < 1e-9){
	  Py(pcnt,tmp[i]) = 1;
	  pcnt++;
	}
      }
    }
    
  } // end ycnt loop
  
  delete[] xlvl;
  delete[] ylvl;
  delete[] tmp;
}


void LagMono2DTransformsCartesian(const int order, const int msh_lin, const fullMatrix<scalar> Px, const fullMatrix<scalar> Py, fullMatrix<scalar> &Lag2MonoX, fullMatrix<scalar> &MonoX2MonoY, fullMatrix<scalar> &MonoY2Lag){
  /* This function returns the 1D transforms for 2D transforms of
     Lagrangian basis. Ax = Lag2MonoX * U, Ay = MonoX2MonoY * Ax,
     U = MonoY2Lag Ay. */

  int N_s = (order+1)*(order+1);
  // 1D basis: phi and integration points
  const polynomialBasis *basis = polynomialBases::find (msh_lin);
  fullMatrix<double> points, weight;
  gaussIntegration::getLine(order*2+1, points, weight);
  int N_G1D = points.size1();
  int N_s1D = order+1;
  fullMatrix<scalar> phi (N_G1D,N_s1D); 
  fullMatrix<double> phiD (N_G1D,N_s1D); 
  basis->f (points, phiD);
  for(int g = 0; g < N_G1D; g++)
    for(int i = 0; i < N_s1D; i++)
      phi(g,i) = (scalar)phiD(g,i);

  // Vandermonde matrix and inverse
  fullMatrix<scalar> V;
  fullMatrix<scalar> Vinv;
  monovandermonde1d(order, points, V);
  V.invert(Vinv);

  // Calculate the complete nodal to modal transform = V1Dinv*phiGL
  fullMatrix<scalar> Transform1D(order+1,order+1);
  Transform1D.gemm(Vinv, phi);
  
  // holds copies of Vinv on the diagonal
  fullMatrix<scalar> DiagVinv(N_s, N_s); 
  for(int p=0;p<order+1;p++)
    for(int i=0;i<order+1;i++)
      for(int j=0;j<order+1;j++)
	DiagVinv(i+p*(order+1),j+p*(order+1)) = Transform1D(i,j);
  
  // Lag2MonoX = diag(monoVinv,order+1) X Px
  Lag2MonoX.resize(N_s,N_s);
  Lag2MonoX.gemm(DiagVinv,Px);
  fullMatrix<scalar> MonoX2Lag(N_s,N_s);
  Lag2MonoX.invert(MonoX2Lag); // This breaks for floats

  // Lag2MonoY = diag(monoVinv,order+1) X Py
  fullMatrix<scalar> Lag2MonoY(N_s,N_s);
  Lag2MonoY.gemm(DiagVinv,Py);
  MonoY2Lag.resize(N_s,N_s);
  Lag2MonoY.invert(MonoY2Lag); // This breaks for floats
    
  // MonoX2MonoY = Lag2MonoY X MonoX2Lag
  MonoX2MonoY.resize(N_s,N_s);
  MonoX2MonoY.gemm(Lag2MonoY,MonoX2Lag);
}
