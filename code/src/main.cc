/*!
  \file main.cc
  \brief Main file: sets up and drives the simulations
  \param[in] argc An integer argument count of the command line arguments
  \param[in] argv An argument vector of the command line arguments
  \return integer 0 upon successful completion
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
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
#include "fullMatrix.h"
#include "polynomialBasis.h"
#include "polynomialsJacobi.h"
#include "Gauss.h"
#include "GmshDefines.h"
#include "simpleMesh.h"
#include "philmesh.h" //PEJ 11/06/2017
#include "scalar_def.h"
#include "dg_functions.h"
#include "deck.h"
#include "init_cond.h"
#include "error_quant.h"
#include "rk.h"
#include "misc.h"
#include "misc_cuda.h"
#include "limiting.h"
#include "dg_solver.h"
#include "communicator.h"
#include "printer.h"
#include "sensor.h"
#include "mem_counter.h"
#include "timers.h"
#include "lagrange_particles.h"
#include "recovery_tools.h"
#include "mixedform.h"
#include <unistd.h> //for sleep function, which causes execution delay (for parallel debug)


//
// Function prototypes
//
void get_element_types(const int order, int &msh_hex, int &msh_qua, int &msh_tri, int &msh_lin){
  if      (order==0)  {msh_hex = MSH_HEX_8; msh_qua = MSH_QUA_4;    msh_tri = MSH_TRI_3;    msh_lin = MSH_LIN_2;  }
  else if (order==1)  {msh_hex = MSH_HEX_8; msh_qua = MSH_QUA_4;    msh_tri = MSH_TRI_3;    msh_lin = MSH_LIN_2;  }
  else if (order==2)  {msh_hex = MSH_HEX_27; msh_qua = MSH_QUA_9;    msh_tri = MSH_TRI_6;    msh_lin = MSH_LIN_3;  }
  else if (order==3)  {msh_hex = MSH_HEX_64; msh_qua = MSH_QUA_16;   msh_tri = MSH_TRI_10;   msh_lin = MSH_LIN_4;  }
  else if (order==4)  {msh_hex = MSH_HEX_125; msh_qua = MSH_QUA_25;   msh_tri = MSH_TRI_15;   msh_lin = MSH_LIN_5;  }
  else if (order==5)  {/*msh_hex = MSH_HEX_216; */msh_qua = MSH_QUA_36;   msh_tri = MSH_TRI_21;   msh_lin = MSH_LIN_6;  }
  else if (order==6)  {/*msh_hex = MSH_HEX_343; */msh_qua = MSH_QUA_49;   msh_tri = MSH_TRI_28;   msh_lin = MSH_LIN_7;  }
  else if (order==7)  {/*msh_hex = MSH_HEX_512; */msh_qua = MSH_QUA_64;   msh_tri = MSH_TRI_36;   msh_lin = MSH_LIN_8;  }
  else if (order==8)  {/*msh_hex = MSH_HEX_729; */msh_qua = MSH_QUA_81;   msh_tri = MSH_TRI_45;   msh_lin = MSH_LIN_9;  }
  else if (order==9)  {/*msh_hex = MSH_HEX_1000; */msh_qua = MSH_QUA_100;  msh_tri = MSH_TRI_55;   msh_lin = MSH_LIN_10; }
  else if (order==10) {/*msh_hex = MSH_HEX_1331; */msh_qua = MSH_QUA_121;  msh_tri = MSH_TRI_66;   msh_lin = MSH_LIN_11; }
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
  //Some diagnostic options:
  int verbose = 0; //if set to 1, there will be A LOT of printing to the screen, helpful for debugging
  int RecoTest = 0; //test recovery procedures and CGR+ICB structures with fake DOF
  int MapTest = 0; //test the interface collocation procedure
  int AuxSolveTest = 0; //test the auxiliary solve for gradient approximation
  int SigFaceTest = 0; //test the auxiliary solve for gradient approximation

  //Method options: It would be good to pass these through run.py eventually,
  /*
    A guide for the unfortunate soul using this code: HAG is a special scheme
    design where I want UhCommon to be the true recovered solution when solving
    for element sigma variable, but I want the interface gradient to be the
    same as what we see in the BR2 method (UhCommon is average instead of recovered).
    So, when HAG==1, I send an averaging operator in to SigFaceMatrices routine
    instead of recovery operator.

    GradMethod: This command affects how the recovery operator itself is formed;
    setting it to 0 builds true recovery operator, setting to 1 builds
    simple averaging operator.

    Rule of thumb: On Cartesian (or really friendly simplex) mesh, set GradMethod_Vis=0.
    Otherwise, set it to 1 to maximize robustness. Chi_Vis=2 is usually appropriate;
    set it higher if the code misbehaves on unstructured grid.

PHIL WAS HERE
   */
  scalar DomainLength = 2.0*acos(-1); //This needs to match domain length; influences IC and also some mesh adjustments for periodic case
  int RKchoice = 3; //order of RK time integration scheme
  int HAG = 0; //So, If HAG==1, then ues avg for interface grad regardless of PSIxR structure
  int GradMethod_Vis = 1; //0 for CGR-R, 1 for BR2, 2 for averaged-biased-recovery
  scalar Chi_Vis = 2.0; //chi parameter to control connectivity in interface gradient solve. 
  
  //Scheme options to solve for gradient when artificial dissipation is applied:
  int GradMethod_AD = GradMethod_Vis; //0 for CGR-R-Naive, 1 for BR2-Naive (naive means that element interior gradient is the broken gradient)
  scalar Chi_AD = Chi_Vis; //jump penalization in interface gradient solve, AD operations

  /*
    Chi_Vis and Chi_AD: the mixedform.cc routine is built for the CGR-R method,
    where the interface correction is Chi*(f-U^-). For BR2 and CGR-J, the correction
    is often written in my work (Chi/2)*(U^+-U^-). Now, in mixedform.cc,
    if PSIxR_global is configured to give the interface average, then the correction
    is instead set to Chi*(Avg-U^-)=(Chi)*(U^+/2 - U^-/2 - U^-)=(Chi/2)*(U^+-U^-).
    So, while it may look strange, the SigFaceMatrices routine in mixedform.cc
    needs no modification to move between CGR-R and BR2
  */

  //Set the epsGen, Beta_S, and Mew_S parameters
  scalar Beta_S = 0.5; //Artificial viscosity strength. 1/2 works well for Cartesian,p1andp3
  scalar Mew_S = 0.1;//0.33; //C spreading strength; 1 for original Reisner apprach; strong effect on timestep stability constraint, because C equation is a high-diffusivity diffusive law.
  scalar epsGen = pow(10,-11); //general "don't divide by zero" epsilon for AD
  scalar Cthresh = pow(10,-8); //if elemCmax>Cthresh, then AD flux is evaluated in, on, and around the element
  //If Cthresh is negative, then AD is always active, independent of sensor

  //code uses quadrature resolution of order*QR_A + QR_B; set the two parameters here
  int QR_A = 2;
  int QR_B = 1;

  //Command to execute the periodcity fix for the HiOCFD5 vortex case.
  //I think both have the power to mess up the simulation outside of that test case.
  int PeriCo2D = 0; 
  int AdjustInterfaceMap = 0; //set to 1 for simplex meshes, HiOCFD5 vortex case.
  
  //END OF PARAMETER POPULATION, I HOPE+++++++++++++++++++++++++++++++

  //////////////////////////////////////////////////////////////////////////   
  //
  // Memory counter init
  //
  //////////////////////////////////////////////////////////////////////////
  MEM_COUNTER mem_counter;
  
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
 
  //////////////////////////////////////////////////////////////////////////   
  //
  // Let's get started!
  //
  //////////////////////////////////////////////////////////////////////////
  if(myid==0){
    printf("\n\n\n");
    printf("Deploy the algorithm...\n");
    printf("\n\n");
    printf("+--------------------------------------------------------------------------------+\n");
    printf("|+------------------------------------------------------------------------------+|\n");
    printf("||                                                                              ||\n");
    printf("||                                                                              ||\n");
    printf("||                                                                              ||\n");
    printf("||                              ALGORITHM DEPLOYED.                             ||\n"); 
    printf("||                                                                              ||\n");
    printf("||                                                                              ||\n");
    printf("||                                                                              ||\n");
    printf("|+------------------------------------------------------------------------------+|\n");
    printf("+--------------------------------------------------------------------------------+\n");
    printf("\n\n\n");
  }
  
  ////////////////////////////////////////////////////////////////////////////
  //
  // Timer stuff
  //
  ////////////////////////////////////////////////////////////////////////////
  TIMERS timers(myid);
  timers.start_timer(0);
  
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
  if(myid==0){printf("Using BLAS\n");}

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
#elif ROP
  if(myid==0){printf("Using ROP (modified Roe) solver, only available for 2d singlefluid\n");}
#elif UPW
  if(myid==0){printf("Using Exact Upwind Flux, only available for scalar advection\n");}
#elif CEN
  if(myid==0){printf("Using Central (average) Flux, only available for scalar advection\n");}
#elif SLAU
  if(myid==0){printf("Using SLAU, only available for singlefluid\n");}
#endif


  //Get the viscosity rule
#ifdef NOVIS
  if (myid == 0){printf("Using NOVIS (viscosity=diffusivity=0)\n");}
#elif CONSTANTVIS
  if (myid == 0){printf("Using CONSTANTVIS (uniform viscosity set in physics/constants.h)\n");}
#elif SUTHERLAND
  if (myid == 0){printf("Using SUTHERLAND (viscosity parameters set in physics/constants.h)\n");}
#endif

  //Get the advective reconstruction approach
#ifdef SIMPLE
  if (myid == 0){printf("Using SIMPLE: typical DG trace values for Riemann solver inputs\n");}
  #elif ICBN
  if (myid == 0){printf("Using ICBN: biased recovery-based reconstructions for Riemann solver inputs\n");}
#endif
  
  // Get the mesh
  std::string fileName = inputs.getMeshfile();
  
  // Setup the limiting
  int limiterMethod = 0;
  if      (inputs.getLimiter() == "hrl")   {limiterMethod = 1; if(myid==0){printf("Using HR limiting (global)\n");}}
  else if (inputs.getLimiter() == "m2l")   {limiterMethod = 2; if(myid==0){printf("Using m2 limiting (global)\n");}}
  else if (inputs.getLimiter() == "hri")   {limiterMethod = 3; if(myid==0){printf("Using HR limiting (per element)\n");}}
  else if (inputs.getLimiter() == "m2i")   {limiterMethod = 4; if(myid==0){printf("Using m2 limiting (per element)\n");}}
  else if (inputs.getLimiter() == "prl")   {limiterMethod = 5; if(myid==0){printf("Using primitive reconstruction limiting (global)\n");}}
  else if (inputs.getLimiter() == "pri")   {limiterMethod = 6; if(myid==0){printf("Using primitive reconstruction limiting (per element)\n");}}
  else if (inputs.getLimiter() == "p0i")   {limiterMethod = 7; if(myid==0){printf("Using p=0 limiting with a sensor\n");}}
  else{limiterMethod = 0; if(myid==0){printf("No limiting\n");}}

  //  return 0;

  // Setup the initial condition type
  bool HiOWvtx = false; //PEJ Edit 01/24/2017
  bool sodphil = false; //PEJ Edit 01/24/2017
  bool shosher = false; //PEJ Edit 10/25/2017
  bool sodcomp = false; //PEJ Edit 10/10/2017
  bool explode = false; //PEJ Edit 10/10/2017
  bool kushjet = false; //PEJ Edit 10/10/2017
  bool normvtx = false; //PEJ Edit 02/09/2017
  bool shckvtx = false; //PEJ Edit 02/23/2017
  bool worsvtx = false; //PEJ Edit 10/17/2017
  bool sinphil = false; //PEJ Edit 05/23/2017
  bool sinphilProject = false; //PEJ Edit 07/05/2017
  bool rhobumpProject = false; //PEJ Edit 08/17/2017
  bool normvtxProject = false; //PEJ Edit 09/08/2017
  bool HiOWvtxProject = false; //PEJ edit 10/02/2017
  bool worsvtxProject = false; //PEJ edit 10/02/2017
  bool explodeProject = false; //PEJ Edit 10/27/2017
  bool tgrvrtx = false; //PEJ edit 09/26/2017 aww yeah
  bool rstrtic = false; //PEJ 10/13/2017: Perform code restart by reading a .pos file
  bool tranvtx = false;
  bool velpert = false;
  bool simplew = false;
  bool sodtube = false;
  bool contact = false;
  bool rhotact = false;
  bool matfrnt = false;
  bool sinegam = false;
  bool expogam = false;
  bool shckint = false;
  bool shuoshe = false;
  bool multint = false;
  bool blast1d = false;
  bool simblst = false;
  bool shckrar = false;
  bool rarecon = false;
  bool sodcirc = false;
  bool rminstb = false;
  bool rmmulti = false;
  bool rtaylor = false;
  bool khdrake = false;
  bool khuramp = false;
  bool khinstb = false;
  bool khblast = false;
  bool khpertu = false;
  bool blastrm = false;
  bool sinephi = false;
  bool sodmono = false;
  bool stffrnt = false;
  bool stfshck = false;
  bool stfbubl = false;
  bool shckdrp = false;
  bool drpwall = false;
  bool jetcrss = false;
  bool prsrflw = false;
  bool injectr = false;
  bool bblwedg = false;
  bool cfplrun = false;
  bool rmawave = false;
  if      (inputs.getInitialCondition()=="simplew") simplew = true;
  else if (inputs.getInitialCondition()=="HiOWvtx") HiOWvtx = true; //PEJ 01/24/2017
  else if (inputs.getInitialCondition()=="sodphil") sodphil = true; //PEJ 01/24/2017
  else if (inputs.getInitialCondition()=="shosher") shosher = true; //PEJ 10/25/2017
  else if (inputs.getInitialCondition()=="sodcomp") sodcomp = true; //PEJ 10/10/2017
  else if (inputs.getInitialCondition()=="explode") explode = true; //PEJ 10/10/2017
  else if (inputs.getInitialCondition()=="kushjet") kushjet = true; //PEJ 10/10/2017
  else if (inputs.getInitialCondition()=="normvtx") normvtx = true; //PEJ 02/08/2017
  else if (inputs.getInitialCondition()=="shckvtx") shckvtx = true; //PEJ 02/23/2017
  else if (inputs.getInitialCondition()=="worsvtx") worsvtx = true; //PEJ 10/17/2017
  else if (inputs.getInitialCondition()=="sinphil") sinphil = true; //PEJ 05/23/2017
  else if (inputs.getInitialCondition()=="sinphilProject") sinphilProject = true; //PEJ 07/05/2017
  else if (inputs.getInitialCondition()=="rhobumpProject") rhobumpProject = true; //PEJ 07/05/2017
  else if (inputs.getInitialCondition()=="normvtxProject") normvtxProject = true; //PEJ 09/08/2017
  else if (inputs.getInitialCondition()=="HiOWvtxProject") HiOWvtxProject = true; //PEJ 10/02/2017
  else if (inputs.getInitialCondition()=="worsvtxProject") worsvtxProject = true; //PEJ 10/02/2017
  else if (inputs.getInitialCondition()=="explodeProject") explodeProject = true; //PEJ 10/27/2017
  else if (inputs.getInitialCondition()=="tgrvrtx") tgrvrtx = true;
  else if (inputs.getInitialCondition()=="rstrtic") rstrtic = true; //PEJ 10/13/2017
  else if (inputs.getInitialCondition()=="tranvtx") tranvtx = true;
  else if (inputs.getInitialCondition()=="velpert") velpert = true;
  else if (inputs.getInitialCondition()=="sodtube") sodtube = true;
  else if (inputs.getInitialCondition()=="contact") contact = true;
  else if (inputs.getInitialCondition()=="rhotact") rhotact = true;
  else if (inputs.getInitialCondition()=="matfrnt") matfrnt = true;
  else if (inputs.getInitialCondition()=="sinegam") sinegam = true;
  else if (inputs.getInitialCondition()=="expogam") expogam = true;
  else if (inputs.getInitialCondition()=="shckint") shckint = true;
  else if (inputs.getInitialCondition()=="shuoshe") shuoshe = true;
  else if (inputs.getInitialCondition()=="multint") multint = true;
  else if (inputs.getInitialCondition()=="blast1d") blast1d = true;
  else if (inputs.getInitialCondition()=="simblst") simblst = true;
  else if (inputs.getInitialCondition()=="shckrar") shckrar = true;
  else if (inputs.getInitialCondition()=="rarecon") rarecon = true;
  else if (inputs.getInitialCondition()=="sodcirc") sodcirc = true;
  else if (inputs.getInitialCondition()=="rminstb") rminstb = true;
  else if (inputs.getInitialCondition()=="rmmulti") rmmulti = true;
  else if (inputs.getInitialCondition()=="rtaylor") rtaylor = true;
  else if (inputs.getInitialCondition()=="khdrake") khdrake = true;
  else if (inputs.getInitialCondition()=="khuramp") khuramp = true;
  else if (inputs.getInitialCondition()=="khinstb") khinstb = true;
  else if (inputs.getInitialCondition()=="khblast") khblast = true;
  else if (inputs.getInitialCondition()=="khpertu") khpertu = true;
  else if (inputs.getInitialCondition()=="blastrm") blastrm = true;
  else if (inputs.getInitialCondition()=="sinephi") sinephi = true;
  else if (inputs.getInitialCondition()=="sodmono") sodmono = true;
  else if (inputs.getInitialCondition()=="stffrnt") stffrnt = true;
  else if (inputs.getInitialCondition()=="stfshck") stfshck = true;
  else if (inputs.getInitialCondition()=="stfbubl") stfbubl = true;
  else if (inputs.getInitialCondition()=="shckdrp") shckdrp = true;
  else if (inputs.getInitialCondition()=="drpwall") drpwall = true;
  else if (inputs.getInitialCondition()=="jetcrss") jetcrss = true;
  else if (inputs.getInitialCondition()=="prsrflw") prsrflw = true;
  else if (inputs.getInitialCondition()=="injectr") injectr = true;
  else if (inputs.getInitialCondition()=="bblwedg") bblwedg = true;
  else if (inputs.getInitialCondition()=="cfplrun") cfplrun = true;
  else if (inputs.getInitialCondition()=="rmawave") rmawave = true;
  else{printf("Invalid initial condition setup. Correct the deck.\n");}

  // Restart option step
  int restart_step = inputs.getRestartStep();
  if((myid==0)&&(restart_step!=0)){printf("Restarting at output step %i.\n",restart_step);}

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
  if (myid == 0) printf("Initial condition declared. Now moving to the mesh loading sequence\n");
  simpleMesh m(myid,numprocs);
  m.load(fileName.c_str());

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
#endif

  


  int elem_type;
  int face_type;
  int msh_hex;
  int msh_qua;
  int msh_tri;
  int msh_lin;
  int nsides; // this offsets j in buildInterfaces function
  int N_N;    // number of neighbors to an element
  scalar refArea; // area of reference element; becomes length in 1D or volume in 3D
  get_element_types(order, msh_hex, msh_qua, msh_tri, msh_lin);
  if     (inputs.getElemType() == "lin"){face_type = MSH_PNT, elem_type = msh_lin; nsides = 0; N_N = 2; refArea = 2.0;}
  else if(inputs.getElemType() == "tri"){face_type = msh_lin, elem_type = msh_tri; nsides = 3; N_N = 3; refArea = 0.5;}
  else if(inputs.getElemType() == "qua"){face_type = msh_lin, elem_type = msh_qua; nsides = 4; N_N = 4; refArea = 4;} 
  else if(inputs.getElemType() == "hex"){face_type = msh_qua, elem_type = msh_hex; nsides = 6; N_N = 6; refArea = 8;}
  else printf("Invalid element type in deck");
  
  
  printf("Set face_type=%d, elem_type=%d, nsides=%d, N_N=%d, refArea=%f\n", face_type, elem_type, nsides, N_N, refArea);

  if (myid == 0) printf("Check3\n"); fflush(stdout);
  if (myid == 0) printf("N_N = %d\n",N_N);
  // Get the nodes, elements, interfaces, normals
  //12/10/2017: Discovered that the entire nodesCONST matrix is stored
  //by each processor; each processor sees the entire mesh,
  //not just the nodes of the elements it owns.
  const fullMatrix<double> &nodesCONST = m.getNodes();
  fullMatrix<double> nodes; nodes.resize(nodesCONST.size1(), nodesCONST.size2());
  for (int i = 0; i < nodes.size1(); i++){
    for (int j = 0; j < nodes.size2(); j++){
      nodes(i,j) = nodesCONST(i,j); }}

  if (myid == 0) printf("Check4\n"); fflush(stdout);
  if (verbose == 1)
    {
      printf("pr=%d: In main: The nodes from getNodes, before periodicity correction (if necessary). nodes(direction, index):\n",myid);
      for (int j = 0; j < nodes.size2(); j++)
	{
	  for (int i = 0; i < nodes.size1(); i++)
	    {
	      //    nodes(i,j) = 0.0;
	      //these are physical node locations read in from the mesh file
	      printf("pr=%d, nodes(%d,%d)=%f,  ",myid, i,j,nodes(i,j));} printf("\n"); }
    }



  nodes = PeriFixNodes(PeriCo2D, nodesCONST); //Make HiOCFD5 vortex mesh adjustment

  //return 0;


   if (verbose == 1 && PeriCo2D>0)
    {
      printf("pr=%d: In main: The nodes from getNodes, after periodicity correction (if necessary). nodes(direction, index):\n",myid);
      for (int j = 0; j < nodes.size2(); j++)
	{
	  for (int i = 0; i < nodes.size1(); i++)
	    {
	      //    nodes(i,j) = 0.0;
	      //these are physical node locations read in from the mesh file
	      printf("pr=%d: nodes(%d,%d)=%f,  ",myid, i,j,nodes(i,j));} printf("\n"); }
    }
   // return 0;
  
   

const std::vector<simpleElement> &elements = m.getElements(elem_type);
 printf("elem_type=%d\n", elem_type);
 if (myid == 0) {printf("Check5\n"); fflush(stdout);}
 /*
   #ifdef USE_MPI
   MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
   printf("\n"); fflush(stdout);
   #endif
 */
 
 
 //  return 0;
 //Get a map of just the N_s master node indices constituting each element:
 std::vector<std::vector<int> > elemNodes;
 elemNodes.resize(elements.size());
 for (int om = 0; om < elements.size(); om++) { elemNodes[om].clear();}
 
 for (int om = 0; om < elements.size(); om++) {
   for (int k = 0; k < elements[om].getNbNodes(); k++){
     elemNodes[om].push_back(elements[om].getNode(k)); }}
 if (verbose > 0)
    {
      printf("pr=%d: The elemNodes structure:\n",myid);
      for (int om = 0; om < elements.size(); om++)
	{
	  printf("pr=%d: om=%d: nodes are {",myid, om);
	  for (int k = 0; k < elements[om].getNbNodes(); k++)
	    {
	      printf("%d, ",elemNodes[om][k]);
	    }
	  printf("}\n");
	}
    }

 // return 0;

  
  nodes = PeriFixNodes_HiOMesh(PeriCo2D, nodes, elemNodes, elem_type, order, N_N); //Make another HiOCFD5 mesh adjustment for p>1 elements
  if (myid == 0) {printf("Returned to main from PeriFixNodes_HiOMesh\n");}

  //m.buildInterfaces(face_type, elem_type,nsides);
  //Build interfaces requires the elements on different partitions
  //to speak to each other, so all partitions
  //need to have their element information properly
  //populated before I call buildInterfaces.
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  m.buildInterfaces(face_type, elem_type, nsides, PeriCo2D); 

  if (myid == 0) {printf("Check6\n"); fflush(stdout);}

 



  //PEJ 11/17/2017: Need an extra routine to check the closure information in interfaces
  //and correct if necessary
  if (AdjustInterfaceMap == 1)
    {
#ifdef USE_MPI
  printf("CATASTROPHE: fixclosures is not programmed for parallel calculations; need to account for partition-boundary interfaces\n");
  MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
  MPI_Finalize();
  return 0;
#endif
      m.fixclosures(face_type, elem_type, nsides, PeriCo2D, DomainLength);
      printf("Check 6.5\n"); fflush(stdout);
    }
  //return 0;
  const std::vector<simpleInterface> &interfaces = m.getInterfaces();
  if (myid == 0) {printf("Check7\n"); fflush(stdout);}
  //m.buildNormals(face_type, elem_type);
 
  m.buildNormals(face_type, elem_type, PeriCo2D, elemNodes, order, N_N);
  if (myid == 0) {printf("Check8\n"); fflush(stdout);}

 

  const fullMatrix<scalar> &normals = m.getNormals();
  if (myid == 0) {printf("Check9\n"); fflush(stdout);}
  m.buildElementMap(elem_type);
  if (myid == 0) {printf("Check10\n"); fflush(stdout);}
  const std::map<int,int> &ElementMap = m.getElementMap();
  if (myid == 0) {printf("Check11\n"); fflush(stdout);}
  /*
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  */
  //Build the ghostElement map and stuff
  m.buildCommunicators(elem_type); // build the indexes to map the ghost elements to my partition
  /*
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  */
  if (myid == 0) {printf("Check12\n"); fflush(stdout);}
  const std::map<int,int> & ghostElementMap = m.getGhostElementMap();
  const std::vector<int> & ghostElement_Key = m.getGhostElement_Key();
  const std::vector<int> & ghostElement_Val = m.getGhostElement_Val();
  const std::vector<int> & ghostElement_HNe = m.getGhostElement_HNe();

  if (verbose > 0)
    {
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
  if (myid == 0) {printf("Check13\n"); fflush(stdout);}
  if (verbose > 0)
    {
      printf("ghostElement_Key/Val/HNe contents for pr=%d: Key is gmsh id of off-partition element, Val is the index where that element's info is stored on pr=%d, and HNe is gmsh id of neighbor of Key/Val element on pr=%d\n", myid, myid, myid);
      for (int j = 0; j < ghostElement_Key.size(); j++)
	{
	  printf("pr=%d, Key=%d, Val=%d, HNe=%d\n", myid, ghostElement_Key[j], ghostElement_Val[j], ghostElement_HNe[j]);
	}
    }
 
  //  if (verbose > 0)
  /*
  {
  printf("ghostElementMap contents for processor %d: first entry is global element id (an element that is off of the processor), second entry is the index where that element's info is stored on pr=%d\n", myid, myid);
  for (std::map<int,int>::const_iterator it=ghostElementMap.begin(); it!=ghostElementMap.end(); ++it)
    {
  printf("pr=%d, ",myid);
  std::cout << it->first << " => " << it->second << '\n';
}
}
  */

  ////////////////////////////////////////////////////////////////////////////   
  //
  // Generer les fonctions de formes, get integration points, weights
  //
  ////////////////////////////////////////////////////////////////////////////   

  // Elements
  const polynomialBasis *basis  = polynomialBases::find(elem_type);  // for the element
  const std::vector<std::vector<int> > &closures = basis->closures;
  fullMatrix<double> points, weight;
  if (verbose > 0){
#ifdef USE_MPI
     MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
#endif
}
  if (myid == 0)
    {
      for (int j = 0; j < closures.size(); j++)
	{
	  printf("basis closure[%d]: nodes are\t{\n",j);
	  for (int m = 0; m < closures[j].size(); m++)
	    {
	      printf("%d, ",closures[j][m]);
	    }
	  printf("}\n");
	}
    }
  
  int QuadRule = order*QR_A + QR_B; //quadrature resolution; order*2+1 is typical but can sometimes be insufficient. 
  if     (inputs.getElemType() == "lin") gaussIntegration::getLine(QuadRule, points, weight);
  else if(inputs.getElemType() == "tri") gaussIntegration::getTriangle(QuadRule, points, weight);
  else if(inputs.getElemType() == "qua") gaussIntegration::getQuad(QuadRule, points, weight);
  else if(inputs.getElemType() == "hex") gaussIntegration::getHexahedron(QuadRule, points, weight);
  

  // Faces
  const polynomialBasis *basisF; // for the edges
#ifdef TWOD
  basisF = polynomialBases::find (msh_lin);
#elif THREED
  basisF = polynomialBases::find (msh_qua);
#endif

  fullMatrix<double> pointsF, weightF;
#ifdef ONED
  pointsF.resize(1,3); weightF.resize(1,1); weightF(0,0)=1;
#elif TWOD
  gaussIntegration::getLine(QuadRule, pointsF, weightF);
#elif THREED
  gaussIntegration::getQuad(QuadRule, pointsF, weightF);
#endif
  
  if (myid == 0 && verbose == 1)
    {
      printf("---Quadrature information:---\n");
      printf("---Element Quadrature:---\n");
      double sumW = 0.0;
      for (int g = 0; g < points.size1(); g++)
	{
	  printf("weight[%d]=%f; ref cords are (",g, weight(g,0));
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",points(g,a));
	    }
	  printf(")\n");
	  sumW += weight(g,0);
	}
      printf("sum of the volume quadrature weights = %f\n", sumW);
      sumW = 0.0;
      printf("---Face Quadrature---:\n");
      for (int g = 0; g < pointsF.size1(); g++)
	{
	  printf("weightF[%d]=%f; ref cords are (",g, weightF(g,0));
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",pointsF(g,a));
	    }
	  printf(")\n");
	  sumW += weightF(g,0);
	}
      printf("sum of the face quadrature weights = %f\n", sumW);
    }
 

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
#elif THREED
  M_s = pow(order+1, 2); //assume hex elements for now, with quad faces
#endif
  int N_T = basis->numFaces;          // number of faces per element            
  int N_E = elements.size();          // number of elements                     (e index)
  int M_T = interfaces.size();        // number of faces                        (t index)
  int N   = N_E*N_s;                  // number of dof for a DG
  int N_G = points.size1();           // number of integration points           (g index)
  int M_G = pointsF.size1();          // number of integration points on a face (g index)
  int N_ghosts = m.getNbGhosts();
  

  /*
#ifdef USE_MPI
     MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
     sleep(1);
#endif
  */
#ifdef USE_MPI
  //placing a barrier here because I want to see the problem setup parameters for all partitions
     MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
#endif

     printf("\nProblem Parameters for processor = %d: elements.size()=%ld\n", myid,  elements.size());
     printf("pr=%d: N_s      = solution points per element              = %i\n",myid, N_s);
     printf("pr=%d: M_s      = supported nodes per face of ref element  = %i\n",myid, M_s);
     printf("pr=%d: N_T      = faces per element                        = %i\n",myid, N_T);
     printf("pr=%d: N_E      = element count                            = %i\n",myid, N_E);
     printf("pr=%d: M_T      = total interface count (with boundaries)  = %i\n",myid, M_T);
     printf("pr=%d: N_G      = quadrature nodes per element             = %i\n",myid, N_G);
     printf("pr=%d: M_G      = quadrature nodes per interface           = %i\n",myid, M_G);
     printf("pr=%d: N_ghosts = count of ghost elements (for partition connections)         = %i\n",myid, N_ghosts);
     printf("pr=%d: D        = number of spatial dimensions             = %d\n",myid, D);
     printf("pr=%d: DF       = interface spatial dimensions             = %d\n", myid, DF);
     if (verbose > 0)
       {
	 printf("This %d processor's element information:\n",myid);
	 for (int om = 0; om < elements.size(); om++)
	   {
	     printf("pr=%d, om=%d\n",myid, om);
	     printf("om=%d: element id is %d\n",om, elements[om].getId());
	     printf("\tphysical tag = %d\n",elements[om].getPhysicalTag());
	     printf("\tpartition = %d\n",elements[om].getPartition());
	     printf("\tnumber of nodes = %d\n", elements[om].getNbNodes());
	     for (int j = 0; j < elements[om].getNbNodes(); j++)
	       {
		 int hold = elements[om].getNode(j);
		 printf("\t\tnode[%d] index = %d, physical location = (%f,%f,%f)\n", j, hold,nodes(0,hold),nodes(1,hold),nodes(2,hold));
	       }
	     printf("---\n");
	   }
       }
       
     //  return 0;

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
#elif THREED
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
  if (myid == 0 && verbose == 1)
    {
      printf("Built the basis on reference element\n");
      for (int k = 0; k < N_s; k++)
	{
	  printf("Basis index %d out of %d:\n",k,N_s);
	  for (int g = 0; g < N_G; g++)
	    {
	      printf("\tquadr. node %d: phi=%f, grad_ref_phi=(",g, phi(g,k));
	      for (int a = 0; a < D; a++)
		{
		  printf("%f, ",dphi(g*D+a, k));
		}
	      printf(")\n");
	    }
	}
      printf("Also built the basis on the reference interface\n");
      for (int k = 0; k < M_s; k++)
	{
	  printf("Basis index %d out of %d:\n",k,M_s);
	  for (int g = 0; g < M_G; g++)
	    {
	      printf("\tquadr. node %d: psi=%f, grad_ref_psi=(",g, psi(g,k));
	      for (int a = 0; a < DF; a++)
		{
		  printf("%f, ",dpsi(g*DF+a, k));
		}
	      printf(")\n");
	    }
	}
    }
  // return 0;
  
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
  //11/06/2017: XYZNodes must adhere to periodicity fix for HiOCFD5 vortex case;
  //no adjustment necessary in XYZnodes, because it uses nodes(,), which I have
  //already corrected..
  //Also, while each processor holds the _nodes array for the
//entire mesh, and node information (like XYZnodes) with an N_E
//in the allocation is only knowledgable of on-partition elements.
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
  //return 0;
  // Element centroids
  //Thus subroutine works with N_N for 1D and 2D but in 3D, count of
  //vertex nodes per element does not match number of faces per element
  int arg2Cen = N_N;
  if (D == 3)
    {
      //3D problem, send an 8 instead, assuming hex elements
      arg2Cen = 8;
    }
  fullMatrix<scalar> XYZCen = m.getElementCentroids(N_E, arg2Cen, XYZNodes);


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

  if (verbose == 1)
    {
      printf("Examining Physical Coordinates, processor number %d:\n", myid);
      printf("---We now have some physical coordinates:---\n");
      printf("---Element coordinates:---\n");
      for (int om = 0; om < N_E; om++)
	{
	  printf("\tpr=%d, om=%d: Centroid= (",myid,om);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",XYZCen(om,a));
	    }
	  printf(")\n");
	  for (int k = 0; k < N_s; k++)
	    {
	      printf("solution node %d: XYZ = (",k);
	      for (int a = 0; a < D; a++)
		{
		  printf("%f, ",XYZNodes(k,om*D+a));
		}
	      printf(")\n");
	    }
	}
      //Coordinates need to match for nodes on
      //opposing sides of the interface; if not, 
      //then through related events, the 
      //interface solutions sent to Riemann solver
      //will be wrong
      printf("---Interface Coordinates---\n");
      for (int t = 0; t < M_T; t++)
	{
	  printf("\tpr=%d, t=%d\n", myid, t);
	  for (int k = 0; k < M_s; k++)
	    {
	      printf("side = 0, solution node %d: XYZ = (",k);
	      for (int a = 0; a < D; a++)
		{
		  printf("%f, ",XYZNodesF(k,(t*2+0)*D+a));
		}
	      printf(")\n");
	    }
	  for (int k = 0; k < M_s; k++)
	    {
	      printf("side = 1, solution node %d: XYZ = (",k);
	      for (int a = 0; a < D; a++)
		{
		  printf("%f, ",XYZNodesF(k,(t*2+1)*D+a));
		}
	      printf(")\n");
	    }
	}
    }
    
  

  //////////////////////////////////////////////////////////////////////////   
  //
  // Build the boundary map (must be done after the normals)
  //
  //////////////////////////////////////////////////////////////////////////
  m.buildBoundary();
  int M_B = m.getBoundarySize();
  printf("pr=%d, M_B %i\n",myid,M_B);
  


  
  
  
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

  m.buildNeighbors(N_N, N_E); //I don't think we need barrier before this because
    
  if (myid == 0 && verbose == 1) {printf("Done with buildNeighbors, returned to main.cc\n"); fflush(stdout);}

 

  //return 0;

  //////////////////////////////////////////////////////////////////////////   
  //
  // Monomial to Lagrange basis transforms
  //
  //////////////////////////////////////////////////////////////////////////
#ifdef ONED
  //PEJ 10/16/2017: monoV has to be inverted, which suggests
  //that it needs to be square; such is not the case 
  //for #points>N_s; so, I need a separate
  //gaussian quadrature set here that matches 
  //N_s for a given order p.
  fullMatrix<double> pointsLIM, weightLIM;
  gaussIntegration::getLine(order*2+1, pointsLIM, weightLIM);
  //The basis"phi" used for DG code as a whole
  //is also not approptiate here, so I will introduce phiLIM
  fullMatrix<scalar> phiLIM (N_G,N_s); 
  fullMatrix<double> phiDLIM (N_G,N_s); 
 
  basis->f (pointsLIM, phiDLIM);
  for(int g = 0; g < pointsLIM.size1(); g++){
    for(int i = 0; i < N_s; i++){
      phiLIM(g,i) = (scalar)phiDLIM(g,i);
    }
  }  
  fullMatrix<scalar> monoV;
  fullMatrix<scalar> monoVinv;
  // monovandermonde1d(order, points, monoV);
  monovandermonde1d(order, pointsLIM, monoV);
  monoV.invert(monoVinv);

  // Go from lagrange to monomial basis (in ref space)
  fullMatrix<scalar> Lag2Mono(N_s,N_s);
  fullMatrix<scalar> Mono2Lag(N_s,N_s);
  //Lag2Mono.gemm(monoVinv, phi);   // Calculate the complete nodal to modal transform = V1Dinv*phiGL
  Lag2Mono.gemm(monoVinv, phiLIM);   // Calculate the complete nodal to modal transform = V1Dinv*phiGL
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
    cartesian_permutations(order, XYZNodes, Px,Py);
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
  //     getPowersXYZG(N_E, N_s, N_G, N_N, M_B, order, XYZG, XYZCen, m.getNeighbors(), shifts, h_powersXYZG);}
  //   else if(inputs.getElemType() == "qua"){
  //     h_powersXYZG = new scalar[L2Msize1*N_G*(N_N+1)*N_E];
  //     getPowersXYZG(N_E, L2Msize1, N_G, N_N, M_B, 2*order, XYZG, XYZCen, m.getNeighbors(), shifts, h_powersXYZG);}
  // }

#elif THREED
  //Don't know how 3D limiting will work yet, but I'd like to declare a few structures 
  //just so the code can run

  printf("WARNING: Lag2Mono, etc. not properly set up in 3D\n");
  fullMatrix<scalar> Lag2MonoX;
  fullMatrix<scalar> MonoX2MonoY;
  fullMatrix<scalar> MonoY2Lag;
  fullMatrix<scalar> monoV;
  /*
  if(cartesian){
    fullMatrix<scalar> Px;
    fullMatrix<scalar> Py;
    cartesian_permutations(order,XYZNodes, Px,Py);
    LagMono2DTransformsCartesian(order, msh_lin, Px, Py, Lag2MonoX, MonoX2MonoY, MonoY2Lag);
    monovandermonde1d(order, pointsF, monoV);
  }
  */
#endif
  if (myid == 0) printf("Passed Marc's Limiting supplements population\n"); fflush(stdout);
  // return 0;
    
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
  // Communication setup. Doing it here because I need XYZNodes_extented
  //
  //////////////////////////////////////////////////////////////////////////
  COMMUNICATOR communicator(N_ghosts, N_s, m, timers, mem_counter);

  //on top of N_E elements, I have the ghost elements.

  //int Ne_AUG = ElementMap.size() + ghostElementMap.size();
  int Ne_AUG = N_E + N_ghosts;
  printf("processor %d: Ne_AUG=%d\n",myid, Ne_AUG);
  //With communicator in place, build the extended XYZNodes array;
  //it holds solution node coordinates for eacl processor's flesh and ghost elements
  scalar* XYZNodes_extended = new scalar[Ne_AUG*D*N_s]; makeZero(XYZNodes_extended,N_s*(N_E+N_ghosts)*D); mem_counter.addToCPUCounter(Ne_AUG*D*N_s*sizeof(scalar));
  //populate XYZNodes_extended for flesh elements of the processor:
  for (int e = 0; e < N_E; e++){
    for (int a = 0; a < D; a++){
      for (int k = 0; k < N_s; k++) {
	XYZNodes_extended[e*D*N_s + a*N_s + k] = XYZNodes(k, e*D + a);
      }
    }
  }
  printf("processor %d passed XYZNodes flesh assignment\n",myid);

  //Use communicator to fill in the coordinates for ghost elements of the processor:
  //Don't need barrier here because the barrier is contained in the communication routine.
  communicator.CommunicateGhosts(D, XYZNodes_extended); //this is like passing a state variable vector U, except the U is node locations and Nfield is D

  if (verbose > 0) {printf("XYZNodes_extended, processor=%d: N_E=%d, Ne_AUG=%d\n", myid, N_E, Ne_AUG);}
  for (int e = 0; e < Ne_AUG; e++)
    {
      scalar sum = 0;
      for (int k = 0; k < N_s; k++)
	{
	  for (int a = 0; a < D; a++)
	    {
	      sum += fabs(XYZNodes_extended[e*D*N_s + a*N_s + k]);
	      if (verbose > 0) {printf("pr=%d: e=%d, k=%d, a=%d, Xa=%f\n", myid, e, k, a, XYZNodes_extended[e*D*N_s + a*N_s + k]);}
	    }
	}
      if (sum < pow(10,-10))
	{
	  printf("Problem: on pr=%d, e=%d (Ne=%d, Ne_AUG=%d), all node coordinates are zero\n",myid,e,N_E, Ne_AUG);
	  exit(1);
	}
      
    }
  

/*
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
*/
  //Use the array XYZNode_extended to define its matrix representation:
  fullMatrix<scalar> XYZNodes_extended_Matrix(N_s, Ne_AUG*D);
  for (int k = 0; k < N_s; k++)
    {
  for (int e = 0; e < Ne_AUG; e++)
    {
  for (int a = 0; a < D; a++)
    {
  XYZNodes_extended_Matrix(k, e*D  + a) = XYZNodes_extended[e*D*N_s + a*N_s + k];
}
}
}

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Waitall(0, &request, &status);
#endif

  //////////////////////////////////////////////////////////////////////////   
  //
  // Calculate the jacobian matrix for each integration point.
  // PEJ 11/30/2017: Doing this for flesh and ghost elements
  //
  //////////////////////////////////////////////////////////////////////////

  // Elements
  printf("pr=%d: Declaring Jac, invJac, etc.\n",myid); fflush(stdout);
  fullMatrix<scalar> Jac(Ne_AUG*D,N_G*D); //the full jacobain matrix at each quadrature point, each element
  fullMatrix<scalar> invJac(Ne_AUG*D,N_G*D); // Inverse Jacobian matrix: dxi/dx
  fullMatrix<scalar> J(Ne_AUG,1);            // determinant of the Jacobian
  fullMatrix<scalar> invJ(Ne_AUG,1);         // determinant of the inverse Jacobian
  // printf("Calling dg_jacobians_elements, pr=%d\n",myid); fflush(stdout);
  dg_jacobians_elements(N_G, N_E, Ne_AUG, XYZNodes, XYZNodes_extended_Matrix, dphi, Jac, invJac, J, invJ);
  // printf("Finished dg_jacobians_elements, pr=%d\n",myid); fflush(stdout);
  // Faces
  fullMatrix<scalar> JacF(M_T*2*D,M_G*DF);
  fullMatrix<scalar> JF(M_T*2,1);            // determinant of the Jacobian
  fullMatrix<scalar> invJF(M_T*2,1);         // determinant of the inverse Jacobian
  // printf("Calling dg_jacobians_face, pr=%d\n",myid); fflush(stdout);
  dg_jacobians_face(M_T, XYZNodesF, dpsi, normals, JacF, JF, invJF);
  
  //return 0;
  //PEJ 10/01/2017: Store product of quadrature weights, detJac 
  //PEJ 10/31/2017: Also store the straight determinant, no multiplication by quadrature weights.
  scalar* detJxW_full = new scalar[Ne_AUG*N_G]; mem_counter.addToCPUCounter(Ne_AUG*N_G*sizeof(scalar));
  scalar* detJ_full = new scalar[Ne_AUG*N_G]; mem_counter.addToCPUCounter(Ne_AUG*N_G*sizeof(scalar));
  for (int e = 0; e < Ne_AUG; e++) {
    for (int g = 0; g < N_G; g++) {
      scalar detLocal = 0.0;
      if (D == 1)
	{
	  detLocal = fabs(Jac(e,g));
	}
      else if (D == 2)
	{
	  detLocal = fabs(Jac(e*D + 0,g*D + 0)*Jac(e*D + 1,g*D + 1) - Jac(e*D + 0,g*D + 1)*Jac(e*D + 1,g*D + 0));
	}
      else if (D == 3)
	{
	  //Use Cramer's rule to calculate determinant
	  scalar m11 = Jac(e*D+0,g*D + 0); scalar m12 = Jac(e*D+0,g*D + 1); scalar m13 = Jac(e*D+0,g*D + 2);
	  scalar m21 = Jac(e*D+1,g*D + 0); scalar m22 = Jac(e*D+1,g*D + 1); scalar m23 = Jac(e*D+1,g*D + 2);
	  scalar m31 = Jac(e*D+2,g*D + 0); scalar m32 = Jac(e*D+2,g*D + 1); scalar m33 = Jac(e*D+2,g*D + 2);
	  detLocal= m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 - m22*m31);
	  detLocal = fabs(detLocal);
	}
      detJ_full[e*N_G + g] = detLocal;
      detJxW_full[e*N_G + g] = detLocal*weight(g,0);
    } //end quadrature node loop
  } //end element loop
  
  if (verbose > 0)
    {
      printf("pr=%d: Element Jacobians, full storage:\n",myid);
      for (int e = 0; e < Ne_AUG; e++)
	{
	  printf("pr=%d: element %d:\n", myid, e);
	  for (int g = 0; g < N_G; g++)
	    {
	      printf("detJ_full(e=%d,g=%d) = %f\n", e, g, detJ_full[e*N_G + g]);
	    }
	}
    }



  //return 0;
  
  /*
  printf("JFace:\n");
  for (int t = 0; t < M_T; t++)
    {
      printf("tGlo=%d: JFace=%f\n",t, JF(t*2,0));
    }
  */
  //return 0;
  //////////////////////////////////////////////////////////////////////////   
  // 
  // Calculate the inverse mass matrices
  //
  //////////////////////////////////////////////////////////////////////////

  scalar* h_Minv = new scalar[N_s*N_s*Ne_AUG]; mem_counter.addToCPUCounter(N_s*N_s*Ne_AUG*sizeof(scalar));
  //  printf("Calling dg_inverse_mass_matrix, pr=%d\n", myid); fflush(stdout);
  dg_inverse_mass_matrix(order, elem_type, inputs.getElemType(), N_s, N_E, Ne_AUG, XYZNodes, XYZNodes_extended_Matrix, QuadRule, detJ_full, h_Minv);
 
  printf("pr=%d: Finished dg_jacobians_elements, dg_inverse_mass_matrix,  and g_jacobians_face\n",myid); fflush(stdout);
 

  ////////////////////////////////////////////////////////////////////////
  //
  //  Some extra geometric tools for initialization and error quantification
  //
  ///////////////////////////////////////////////////////////////////////

  //get high-resolution quadrature parameters
  fullMatrix<double> SuperPoints, SuperWeight;
  if     (inputs.getElemType() == "lin") gaussIntegration::getLine(2*QuadRule, SuperPoints, SuperWeight);
  else if(inputs.getElemType() == "tri") gaussIntegration::getTriangle(2*QuadRule, SuperPoints, SuperWeight);
  else if(inputs.getElemType() == "qua") gaussIntegration::getQuad(2*QuadRule, SuperPoints, SuperWeight);
  else if(inputs.getElemType() == "hex") gaussIntegration::getHexahedron(2*QuadRule, SuperPoints, SuperWeight);

  int GQsuper = SuperPoints.size1();           // number of Super-integration points   (g index)
  if (verbose > 0) {printf("GQsuper=%d, pr=%d\n", GQsuper,myid);}
  //Get the solution/testing basis on reference element, super-quadrature points
  fullMatrix<scalar> SuperPhi (GQsuper,N_s); 
  fullMatrix<double> SuperPhiD (GQsuper,N_s); 
  basis->f (SuperPoints, SuperPhiD);
  for(int g = 0; g < GQsuper; g++){
    for(int i = 0; i < N_s; i++){
      SuperPhi(g,i) = (scalar)SuperPhiD(g,i);
    }
  } 
  //Also, get ref element solution basis gradient, super-quadrature, for sake of Jac calculation
  fullMatrix<scalar> SuperdPhi (GQsuper*D, N_s);
  for(int g = 0; g < GQsuper; g++){
    basis->df(SuperPoints(g,0),SuperPoints(g,1),SuperPoints(g,2),grads);
    for(int alpha = 0; alpha < D; alpha ++){
      for(int i = 0; i < N_s; i++){
  	SuperdPhi(g*D+alpha,i) = (scalar)grads[i][alpha];  //see paper for indexing p.6
      }	  
    }    
  }

  //Get the physical coordinates of the super-resolution quadrature points
  fullMatrix<scalar> SuperXYZ_GQ (GQsuper, N_E*D);
  SuperXYZ_GQ.gemm( SuperPhi, XYZNodes);

  //One more thing: the high-resolution jacobian. I'm going to get it for 
  //every quadrature point in preparation for future exploits.
  //Copying some stff from dg_jacobians_elements in dg/dg_functions.cc.
  //Need to do this in subroutine eventually
  fullMatrix<scalar> SuperdetJ (N_E, GQsuper);
  fullMatrix<scalar> SuperJac(N_E*D,GQsuper*D);
// fullMatrix<scalar> SuperinvJac(N_E*D,GQsuper*D); //why not? Because that's a lot of memory
  SuperJac.gemm(XYZNodes.transpose(),SuperdPhi.transpose());
  scalar det_local = 0.0;
  for (int e = 0; e < N_E; e++)
    {
      for (int g = 0; g < GQsuper; g++)
	{
#ifdef ONED
	  det_local = sqrt(SuperJac(e,g)*SuperJac(e,g));
#endif
#ifdef TWOD
	  //not sure about this indexing yet,07/05
	  det_local = SuperJac(e*D+0,g*D+0)*SuperJac(e*D+1,g*D+1)-SuperJac(e*D+1,g*D+0)*SuperJac(e*D+0,g*D+1);
#endif
#ifdef THREED
	  //Use cramer's rule to get the determinant
	  {
	    scalar m11 = SuperJac(e*D+0,0); scalar m12 = SuperJac(e*D+0,1); scalar m13 = SuperJac(e*D+0,2);
	    scalar m21 = SuperJac(e*D+1,0); scalar m22 = SuperJac(e*D+1,1); scalar m23 = SuperJac(e*D+1,2);
	    scalar m31 = SuperJac(e*D+2,0); scalar m32 = SuperJac(e*D+2,1); scalar m33 = SuperJac(e*D+2,2);
	    det_local = m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 - m22*m31);
	  }
#endif
	  SuperdetJ(e,g) = fabs(det_local);
	  //Can calculate the SuperInvJac as well,
	  //but I prefer not to because I have no need of it
	  //and the memory cost is severe.
	  /*
#ifdef ONED
	  SuperinvJac(e,g) =  1.0/SuperJac(e,g);
#endif
#ifdef TWOD
	  SuperinvJac(e*D+0,g*D+0) =  1.0/det_local*SuperJac(e*D+1,g*D+1);
	  SuperinvJac(e*D+1,g*D+0) = -1.0/det_local*SuperJac(e*D+1,g*D+0);
	  SuperinvJac(e*D+0,g*D+1) = -1.0/det_local*SuperJac(e*D+0,g*D+1);
	  SuperinvJac(e*D+1,g*D+1) =  1.0/det_local*SuperJac(e*D+0,g*D+0);
#endif
#ifdef THREED
	  fullMatrix<double> JacNode(D,D);
	  fullMatrix<double> invJacNode(D,D);
	  for (int a1 = 0; a1 < D; a1++)
	    {
	  for (int a2 = 0; a2 < D; a2++)
	    {
	  //Not sure if this indexing is right
	  JacNode(a1,a2) = SuperJac(e*D + a1, g*D + a2);
	}
	}
	  //Invert JacNode, store result in invJacNode
	  JacNode.invert(invJacNode);
	  //Now, relay the result to global storage,
	  //using same jumbled ordering as 2D case
	  for (int a1 = 0; a1 < D; a1++)
	    {
	  for (int a2 = 0; a2 < D; a2++)
	    {
	  //not sure if this indexing is right
	  SuperinvJac(e*D + a2, g*D + a1) = invJacNode(a1,a2);
	}
	}
#endif
*/
	}
    }
  printf("pr=%d: Done building SuperJac and SuperInvJac for initialization and error quantification. Next up, initialize the field variables.\n",myid);
  //return 0;
  //////////////////////////////////////////////////////////////////////////   
  //
  // Initialize the unknowns
  //
  //////////////////////////////////////////////////////////////////////////
/*
#ifdef USE_MPI
  //Another waitstation
  MPI_Barrier(MPI_COMM_WORLD);
  //  MPI_Waitall(0,&request, &status);
#endif
*/
  fullMatrix<scalar> U(N_s, N_E*N_F);
  fullMatrix<scalar> Us(N_s, N_E*N_F);
  fullMatrix<scalar> Ustar(N_s, N_E*N_F);
#ifdef SCALARAD
  if(sinphil) init_dg_sinphil_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 05/23/2017
  if(sinphilProject) init_dg_sinphilProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, U); //PEJ 07/05/2017
#endif
#ifdef SINGLEFLUID
  if(tranvtx) init_dg_tranvtx_singlefluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if(HiOWvtx) init_dg_HiOWvtx_singlefluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs()); //PEJ 01/24/2017
  else if(sodphil) 
    {
      printf("Initializing phil's SOD problem\n");
      init_dg_sodphil_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 01/24/2017
    }
  else if(shosher) 
    {
      printf("Initializing Shu-Osher problem, Singlefluid\n");
      init_dg_shosher_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 10/25/2017
    }
  else if(sodcomp) 
    {
      printf("Initializing the compact SOD problem (3 waves left, 3 waves right)\n");
      init_dg_sodcomp_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 01/24/2017
    }
   else if(explode) 
    {
      printf("Initializing the exploding SOD problem \n");
      init_dg_explode_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 01/24/2017
    }
  else if(kushjet) 
    {
      printf("Initializing the ambient atmosphere for Kushner Jet problem \n");
      init_dg_kushjet_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 01/24/2017 
    }
  else if(normvtx)
    {
      printf("Initializing phil's friendly vortex problem\n");
      init_dg_normvtx_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 02/09/2017
    }
  else if (shckvtx)
    {
      printf("Initializing Rault's shock-vortex interaction\n");
      init_dg_shckvtx_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 02/23/2017
    }
  else if (worsvtx)
    {
      printf("Initializing Rault's shock-vortex interaction, Workshop version\n");
      init_dg_worsvtx_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 10/17/2017
    }
  else if (sinphil)
    {
      printf("Initializing phil's sine wave, pr=%d\n",myid);
      init_dg_sinphil_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 05/23/2017
    }
  else if (sinphilProject)
    {
      printf("Initizalizing phil's sine wave through Galerkin projection\n");
      init_dg_sinphilProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, U); //PEJ 07/05/2017
    }
  else if (rhobumpProject)
    {
      printf("Initizalizing phil's rho bump through Galerkin projection\n");
      init_dg_rhobumpProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, U); //PEJ 08/17/2017
    }
  else if (normvtxProject)
    {
      printf("Initizalizing normalized isentropic vortex through Galerkin projection\n");
      init_dg_normvtxProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, U); //PEJ 09/08/2017
    }
  else if (HiOWvtxProject)
    {
      printf("Initializing HiOCFD5 vortex problem through Galerkin projection\n");
      init_dg_HiOWvtxProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, U); //PEJ 10/02/2017
    }
  else if (worsvtxProject)
    {
      printf("Initializing Rault's shock-vortex interaction, Workshop version, through Galerkin projection\n");
      init_dg_worsvtxProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, U); //PEJ 10/17/2017
    }
  else if (explodeProject)
    {
      printf("Initializing square SOD explosion problem, through Galerkin projection\n");
      init_dg_explodeProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, U); //PEJ 10/17/2017
    }
  else if (tgrvrtx)
    {
      printf("Initializing Taylor-Green vortex through solution point equivalence\n");
      init_dg_tgrvrtx_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 05/23/2017
    }
  else if (rstrtic)
    {
      int KUSH_START = 20;
      scalar rsh = 1.4;
      printf("Initializing the DOF by reading in a .pos file, looking in current folder for something that ends in %d.pos\n",KUSH_START);
      printf("Using ratio of specific heats gamma=%f\n",rsh);
      std::vector<std::string> _fnames;
      std::vector<std::string> _names;
      _names.push_back("Rho");   _fnames.push_back("rho");
      _names.push_back("Ux");    _fnames.push_back("ux");
#ifdef TWOD
      _names.push_back("Uy");    _fnames.push_back("uy"); 
#endif
#ifdef THREED
      _names.push_back("Uy");    _fnames.push_back("uy"); 
      _names.push_back("Uz");    _fnames.push_back("uz");
#endif 
      _names.push_back("P");     _fnames.push_back("p");
 
      scalar* solution_start = new scalar[N_E*N_F*N_s];
      //scalar &time_dummy;
      //*time_dummy = 0.0;
      m.readSolutionNoTime(N_s, N_E, elem_type, _fnames, _names, KUSH_START, solution_start);
      for (int e = 0; e < N_E; e++)
	{
	  for (int i = 0; i < N_s; i++)
	    {
#ifdef ONED
	      scalar rho = solution_start[(e*N_F+0)*N_s+i];
	      scalar u = solution_start[(e*N_F+1)*N_s+i];
	      scalar p = solution_start[(e*N_F+(N_F-1))*N_s+i];
	      U(i, e*N_F+0) = rho;
	      U(i, e*N_F+1) = rho*u;
	      U(i, e*N_F+2) = p/(rsh-1.0) + 0.5*rho*(u*u);
#endif
#ifdef TWOD
	      scalar rho = solution_start[(e*N_F+0)*N_s+i];
	      scalar u = solution_start[(e*N_F+1)*N_s+i];
	      scalar v = solution_start[(e*N_F+2)*N_s+i];
	      scalar p = solution_start[(e*N_F+3)*N_s+i];
	      U(i, e*N_F+0) = rho;
	      U(i, e*N_F+1) = rho*u;
	      U(i, e*N_F+2) = rho*v;
	      U(i, e*N_F+3) = p/(rsh-1.0) + 0.5*rho*(u*u + v*v);
#endif
#ifdef THREED
	      scalar rho = solution_start[(e*N_F+0)*N_s+i];
	      scalar u = solution_start[(e*N_F+1)*N_s+i];
	      scalar v = solution_start[(e*N_F+2)*N_s+i];
	      scalar w = solution_start[(e*N_F+3)*N_s+i];
	      scalar p = solution_start[(e*N_F+4)*N_s+i];
	      U(i, e*N_F+0) = rho;
	      U(i, e*N_F+1) = rho*u;
	      U(i, e*N_F+2) = rho*v;
	      U(i, e*N_F+3) = rho*v;
	      U(i, e*N_F+4) = p/(rsh-1.0) + 0.5*rho*(u*u + v*v + w*w);
#endif
	    
	    
	    }
	}
      
      delete[] solution_start;
    }
#elif RADSINGLEFLUID
  if(sodphil) 
    {
      printf("Initializing phil's SOD problem, RAD Singlefluid\n");
      init_dg_sodphil_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 01/24/2017
    }
  else if(shosher) 
    {
      printf("Initializing Shu-Osher problem, RadSinglefluid\n");
      init_dg_shosher_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 10/25/2017
    }
  else if (explode)
    {
      printf("Initializing the explosion problem, RAD Singlefluid\n");
      init_dg_explode_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 01/24/2017
    }
  else if (worsvtx)
    {
      printf("Initializing Rault's shock-vortex interaction, Workshop version\n");
      init_dg_worsvtx_singlefluid(N_s, N_E, XYZNodes, U); //PEJ 10/17/2017
    }
  else if (worsvtxProject)
    {
      printf("Initializing Rault's shock-vortex interaction, Workshop version, through Galerkin projection\n");
      init_dg_worsvtxProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, U); //PEJ 10/17/2017
    }
  else if (explodeProject)
    {
      printf("Initializing square SOD explosion problem, through Galerkin projection\n");
      init_dg_explodeProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, U); //PEJ 10/17/2017
    }
#elif MULTIFLUID
  if     (simplew) init_dg_simplew_multifluid(N_s, N_E, XYZNodes, U);
  else if(velpert) init_dg_velpert_multifluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if(tranvtx) init_dg_tranvtx_multifluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if(sodtube) init_dg_sodtube_multifluid(N_s, N_E, XYZNodes, U);
  else if(sodmono) init_dg_sodmono_multifluid(N_s, N_E, XYZNodes, U);
  else if(contact) init_dg_contact_multifluid(N_s, N_E, XYZNodes, U);
  else if(rhotact) init_dg_rhotact_multifluid(N_s, N_E, XYZNodes, U);
  else if(matfrnt) init_dg_matfrnt_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(sinegam) init_dg_sinegam_multifluid(N_s, N_E, XYZNodes, U);
  else if(expogam) init_dg_expogam_multifluid(N_s, N_E, XYZNodes, U);
  else if(shckint) init_dg_shckint_multifluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getFinalTime());
  else if(shuoshe) init_dg_shuoshe_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(multint) init_dg_multint_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(blast1d) init_dg_blast1d_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(simblst) init_dg_simblst_multifluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if(shckrar) init_dg_shckrar_multifluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if(rarecon) init_dg_rarecon_multifluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if(sodcirc) init_dg_sodcirc_multifluid(N_s, N_E, XYZNodes, U);
  else if(rminstb) init_dg_rminstb_multifluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if(rmmulti) init_dg_rmmulti_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(rtaylor) init_dg_rtaylor_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(khdrake) init_dg_khdrake_multifluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if(khuramp) init_dg_khuramp_multifluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if(khinstb) init_dg_khinstb_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(khblast) init_dg_khblast_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(khpertu) init_dg_khpertu_multifluid(N_s, N_E, XYZNodes, XYZCen, U);
  else if(blastrm) init_dg_blastrm_multifluid(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
#elif PASSIVE
  if (sinephi) init_dg_sinephi_passive(N_s, N_E, XYZNodes, U);
  else if (sodmono) init_dg_sodmono_passive(N_s, N_E, XYZNodes, U);
#elif STIFFENED
  if (stffrnt) init_dg_stffrnt_stiffened(N_s, N_E, XYZNodes, XYZCen, U);
  else if (stfshck) init_dg_stfshck_stiffened(N_s, N_E, XYZNodes, XYZCen, U);
  else if (stfbubl) init_dg_stfbubl_stiffened(N_s, N_E, XYZNodes, XYZCen, U);
  else if (shckdrp) init_dg_shckdrp_stiffened(N_s, N_E, XYZNodes, XYZCen, U);
  else if (drpwall) init_dg_drpwall_stiffened(N_s, N_E, XYZNodes, XYZCen, U);
  else if (jetcrss) init_dg_jetcrss_stiffened(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if (prsrflw) init_dg_prsrflw_stiffened(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if (injectr) init_dg_injectr_stiffened(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if (bblwedg) init_dg_bblwedg_stiffened(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
  else if (cfplrun) init_dg_cfplrun_stiffened(N_s, N_E, XYZNodes, XYZCen, U, m, elem_type, inputs.getInitialConditionInputs());
  else if (rmawave) init_dg_rmawave_stiffened(N_s, N_E, XYZNodes, XYZCen, U, inputs.getInitialConditionInputs());
#endif

  if (order0) average_cell_p0(N_s, N_E, U);
  /*
  printf("Flow initialization complete, pr=%d\n",myid);
  printf("nodal IC, field variable 0:\n");
      for (int e = 0; e < N_E; e++)
	{
	  printf("element %d:\n", e);
	  for (int k = 0; k < N_s; k++)
	    {
	      printf("\t\t@X=(%f,%f,%f): U = %f\n",XYZNodes(k, e*D+0),XYZNodes(k,e*D+1),XYZNodes(k,e*D+2),U(k, e*N_F + 0));
	    }
	    }*/
      // return 0;

  //////////////////////////////////////////////////////////////////////////   
  // 
  // Build the map
  //
  //////////////////////////////////////////////////////////////////////////
  int* h_map = new int[M_s*M_T*N_F*2]; mem_counter.addToCPUCounter(M_s*M_T*N_F*2*sizeof(int));
  int* h_invmap = new int[M_s*N_N*N_E*N_F*2]; mem_counter.addToCPUCounter(M_s*N_N*N_E*N_F*2*sizeof(int));
   //In the initial FaceFromElem, the global interface addresses
  //were stored in order corresponding to sig value of the element.
  //Instead, Alt approach stores the interfaces in the order
  //they make an appearance in invmap
  //05/25/2017: Altered dg_mappings to also populate Alt_FaceFromElem, which tells each element who its interfaces are.
  //Had to spend some time getting the storage order correct, see BuildSigMatrices in mixedform.cc for an adventure
  int* Alt_FaceFromElem = new int[N_E*N_N]; mem_counter.addToCPUCounter(N_E*N_N*sizeof(int));
  // int* Alt_FaceFromElem = new int[(N_E) * (D+N_N)]; mem_counter.addToCPUCounter((N_E)*(D+N_N)*sizeof(int));
  //for (int j = 0; j < N_E*(D+N_N); j++){Alt_FaceFromElem[j] = -1;}

  dg_mappings(myid, M_s, M_T, N_s, N_E, N_N, interfaces, ElementMap, ghostElement_Key, ghostElement_Val, ghostElement_HNe, closures, h_map, h_invmap, Alt_FaceFromElem);
  

 	


 
  //Build a better inverse mapping, where the physical locations of partner
  //nodes on each element of a given interface always match.
  int* Better_InvMap = new int[M_T*N_F*2*M_s]; mem_counter.addToCPUCounter(M_T*N_F*2*M_s*sizeof(int));
  for (int t = 0; t < M_T; t++)
    {
      //Use map[] to identify the two partner nodes
      for (int j = 0; j < M_s; j++)
	{
	  for (int fc =0; fc < N_F; fc++)
	    {
	      //	  int fc=0;
	      int face;
	      for (int d = 0; d < 2; d++) //which side of the interface
		{
		  //face = t*(N_F*2*M_s) + f*(2*M_s) + d*M_s + j;
		  face = ((t*N_F+fc)*2+d)*M_s+j;
		  int serial_element_node_location = h_map[face];
		  /*
		    element_node is the fc=whatever address of the (element,node index)
		    entry corresponding to thig global interface's jth node,
		    but in global storage system.
		  */
	      // Get the interface we will operate on and relate it to the index of U
	      Better_InvMap[face] = serial_element_node_location;
	      //printf("t=%d, Better_InvMap[%d]=%d\n", t, face, serial_element_node_location);
	    }
	}
	}
    }
      
 
  
  //return 0;


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

 

  //PEJ EDIT: Create Space for the pointer format of phigF, other BR2 necessities;
  scalar* h_phigF        = new scalar[M_G*D*N_N*N_s];     makeZero(h_phigF,M_G*D*N_N*N_s);            mem_counter.addToCPUCounter(M_G*D*N_N*N_s*sizeof(scalar));
  scalar* h_dphigF       = new scalar[D*M_G*D*N_N*N_s];   makeZero(h_dphigF,D*M_G*D*N_N*N_s);         mem_counter.addToCPUCounter(D*M_G*D*N_N*N_s*sizeof(scalar));
  scalar* h_invJac_trace = new scalar[M_G*D*D*N_N*N_E*D]; makeZero(h_invJac_trace,M_G*D*D*N_N*N_E*D); mem_counter.addToCPUCounter(M_G*D*D*N_N*N_E*D*sizeof(scalar));
  scalar* h_Jac_trace    = new scalar[M_G*D*D*N_N*N_E*D]; makeZero(h_Jac_trace,M_G*D*D*N_N*N_E*D); mem_counter.addToCPUCounter(M_G*D*D*N_N*N_E*D*sizeof(scalar));
  scalar* h_weight_trace = new scalar[D*N_N*M_G]        ; makeZero(h_weight_trace,D*N_N*M_G)        ; mem_counter.addToCPUCounter(D*N_N*M_G*sizeof(scalar));
  //My BR2 edges/elements/face map: The makeZero argument is apparently not applicable to int
  int* BR2_Map = new int[M_T*4]; mem_counter.addToCPUCounter(M_T*4*sizeof(int));
  
  //Poplulate the BR2_Map array, and FaceFromElem
  //For a given element, FaceFromElem holds the global addresses
  //of the N_N bordering interfaces
  //Also, when a global interface is double-referenced, 
  //FaceFrom Elem needs to hold the lower of the 2 addresses
  int* FaceFromElem = new int[Ne_AUG * N_N];
  for (int j = 0; j < N_E*N_N; j++)
    {
      //once this is NOT -9000, we know that FaceFrom Elem
      //has been populated, and we can leave it alone
      FaceFromElem[j] = -9000;
    }
  
  for (int t = 0; t < M_T; t++) 
    {
      if (verbose > 0) {printf("t=%d, pr=%d\n", t, myid);}
      const simpleInterface &face = interfaces[t]; //fetches the interface corresponding to index t
      const simpleElement *el0 = face.getElement(0);
      const simpleElement *el1 = face.getElement(1);
      if ( ElementMap.find(el1->getId()) == ElementMap.end() ) {
	// not found: the interface's el1 is not in ElementMap
	if (verbose > 0) {printf("For this interface el1 is ghost: el0id=%d, el1id=%d, cl0=%d, cl1=%d, physical=%d\n", ElementMap.at(el0->getId()), ghostElementMap.at(el1->getId()), face.getClosureId(0), face.getClosureId(1), face.getPhysicalTag());}
	BR2_Map[t*4 + 0] = ElementMap.at(el0->getId());
	BR2_Map[t*4 + 2] = m.SeekGhost(el0->getId(), el1->getId());
	//BR2_Map[t*4 + 2] = ghostElementMap.at(el1->getId());
	BR2_Map[t*4 + 1] = face.getClosureId(0);
	BR2_Map[t*4 + 3] = face.getClosureId(1);
	//Populate XYZNodes_extended for recovery procedures.
	//This is fortunately a two-way exchange: both processors
	//need something from the other processor
	/*
	  for (int pr = 0; pr < numprocs; pr++)
	  {
	  if (myid == pr)
	  {
		
	  }
	  }
	  for (int k = 0; k < N_s; k++){
	  for (int a = 0; a < D; a++){
	  XYZNodes_extended(k, BR2_Map[t*4+2]*D+a) = XYZNodes
	*/
      } 
      else {
	// found
	if (verbose > 0) {printf("For this interface el1 is on local partition: el0id=%d, el1id=%d, cl0=%d, cl1=%d, physical=%d\n", ElementMap.at(el0->getId()), ElementMap.at(el1->getId()), face.getClosureId(0), face.getClosureId(1), face.getPhysicalTag());}
	BR2_Map[t*4 + 0] = ElementMap.at(el0->getId());
	BR2_Map[t*4 + 2] = ElementMap.at(el1->getId());
	BR2_Map[t*4 + 1] = face.getClosureId(0);
	BR2_Map[t*4 + 3] = face.getClosureId(1);

	
      }
      //sig_index takes multiple values corresponding to sam face
      int omA = BR2_Map[t*4 + 0];
      int omB = BR2_Map[t*4 + 2];
      int sigA = BR2_Map[t*4 + 1];
      int sigB = BR2_Map[t*4 + 3];
      if (FaceFromElem[omA*N_N + sigA%N_N] == -9000)
	{
	  FaceFromElem[omA*N_N + sigA%N_N] = t;
	}
      if (FaceFromElem[omB*N_N + sigB%N_N] == -9000)
	{
	  FaceFromElem[omB*N_N + sigB%N_N] = t;
	}
      
      //if in parallel and el1=NULL, it means that the interface's other
      //element is in a neighboring partition. Use ghostElementMap instead
      /*
	if (el1 == NULL && numprocs > 1)
	{
	//e11 = -1;
	printf("EL1==NULL: For this interface: el0id=%d, el1id=%d, cl0=%d, cl1=%d, physical=%d\n", ElementMap.at(el0->getId()), -1, face.getClosureId(0), face.getClosureId(1), face.getPhysicalTag());
	}
	else
	{
	printf("For this interface: el0id=%d, el1id=%d, cl0=%d, cl1=%d, physical=%d\n", ElementMap.at(el0->getId()), ElementMap.at(el1->getId()), face.getClosureId(0), face.getClosureId(1), face.getPhysicalTag());
	}
      */
      //printf("For this interface: el0id=%d, el1id=%d, cl0=%d, cl1=%d, physical=%d\n", ElementMap.at(el0->getId()), -1, face.getClosureId(0), face.getClosureId(1), face.getPhysicalTag());
      /*
	for(int d = 0; d < 2; d++)
	{
	const simpleElement *el = face.getElement(d); //grabs the element on a specific side of the interface
	int id = face.getClosureId(d); //ClosureID merits some explanation
	//Quad: Closure ID can be {0..7}: it tells us which face of the reference element we are on; {0..3} run in one direction, {4..7} run in the other.
	//Triangle: Closure ID runs {0..5} instead, because a triangle has only 3 sides
	int el_index = ElementMap.at(el->getId());
	int sig_index = id;// % N_N;
	BR2_Map[t*4 + 2*d + 0] = el_index;
	BR2_Map[t*4 + 2*d + 1] = sig_index;
	printf("pr=%d, t=%d, side=%d, el_index=%d, sig_index=%d\n", myid, t, d, el_index, sig_index);
	  
	//sig_index takes multiple values corresponding to sam face
	if (FaceFromElem[el_index*N_N + sig_index%N_N] == -9000)
	{
	FaceFromElem[el_index*N_N + sig_index%N_N] = t;
	}
	  
	}
      */
      /*
	for (int a = 0; a < D; a++)
	{
	printf("iterface normal(face=%d,dir=%d) = %f\n", t, a, normals(a, t));
	}
      */ 
    }

//Have a look at BR2_Map for each partition
  if (verbose > 0) 
    {
      printf("BR2_Map contents for pr=%d:\n",myid);
      for (int t = 0; t < M_T; t++)
	{
	  printf("processor %d, interface %d: omA=%d, omB=%d, sigA=%d, sigB=%d\n",myid,t,BR2_Map[t*4 + 2*0 + 0],BR2_Map[t*4 + 2*1 + 0],BR2_Map[t*4 + 2*0 + 1],BR2_Map[t*4 + 2*1 + 1]);
	}
    }

  //  return 0;
      /*
    BR2_Map seems to be properly populated and accounts
    for cross-partition traffic. Next Up, need to make
    sure my approach works well through recovery mixed-form approach;
    recovery can wait a bit.
   */
  

  if (verbose == 1)
    {
      printf("in Main: points and weights from first call:\n");
      for (int g = 0; g < points.size1(); g++)
	{
	  printf("weight[%d] = %f\n", g, weight(g,0));
	}
      printf("Quadrature nodes, ref element:\n");
      for (int g = 0; g < points.size1(); g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",points(g,a));
	    }
	  printf(")\n");
	}

      //Phil Edit 05/19/2017
      //Inspect the node locations
      
      for (int e = 0; e < N_E; e++)
	{
	  for (int k = 0; k < N_s; k++)
	    {
	      printf("element %d, node %d = (",e,k);
	      for (int a = 0; a < D; a++)
		{
		  printf("%f,",XYZNodes(k, e*D + a));
		}
	      printf(")\n");
	    }
	  printf("\n");
	}

      //Inspect initial condition:
#ifdef THREED
      printf("nodal IC, field variable 0:\n");
      for (int e = 0; e < N_E; e++)
	{
	  printf("element %d:\n", e);
	  for (int k = 0; k < N_s; k++)
	    {
	      printf("\t\t@X=(%f,%f,%f): U = %f\n",XYZNodes(k, e*D+0),XYZNodes(k,e*D+1),XYZNodes(k,e*D+2),U(k, e*N_F + 0));
	    }
	}
#endif
      printf("invmap:\n");
      for (int e = 0; e < N_E; e++)
	{
	  for (int f = 0; f < 1; f++)
	    {
	      for (int s = 0; s < N_N; s++)
		{
		  for (int m = 0; m < M_s; m++)
		    {
		      printf("invmap[((%d*N_F+%d)*M_s*N_N + %d*M_s + %d)*2 + 0] = %d\n",e,f,s,m,h_invmap[((e*N_F+0)*M_s*N_N +  s*M_s + m)*2+0]);
		    }
		}
	    }
	  printf("\n");
	}

       printf("Better_InvMap:\n");
  for (int t = 0; t < M_T; t++)
    {
      printf("Interface %d:\n",t);
      for (int j = 0; j < M_s; j++)
	{
	  for (int fc = 0; fc < 1; fc++)
	    {
	      printf("\tlocal node %d: partner integers are (",j);
	      for (int d = 0; d < 2; d++)
		{
		  printf("%d ,",Better_InvMap[((t*N_F+fc)*2+d)*M_s+j] );
		}
	      printf(")");
	    }
	  printf("\n");
	}
    }
  

    } //end of verbose == 1 case
  //return 0;
  //Have a look at the BR2_Map array
  if (verbose > 0)
    {
      printf("-----The BR2_Map array:-----\n");
      for (int t = 0; t < M_T; t++)
	{
	  printf("Interface %2d: om_A=%2d, sig_A=%d, om_B=%2d, sig_B=%d\n",t,BR2_Map[t*4+0],BR2_Map[t*4+1],BR2_Map[t*4+2],BR2_Map[t*4+3]);
	}
    }



   //PEJ 05/19/2017: Use element and iterface information to build Recovery PSIxR matrix at each interface.
  //08/17/2017: Also building the biased PSIxR operators in here.
  //10/18/2017: PSIxR_global and PSIxR_global_AD can apply different approaches to get interface state
  //11/27/2017: For parallel case, inverface-based recovery works on >N_E elements per partition;
 
  
  int N_reco = 2*N_s; //dimension of recovery basis 
  scalar* PSIxR_global = new scalar[M_T*M_G*N_reco]; mem_counter.addToCPUCounter(M_T*M_G*N_reco*sizeof(scalar));

  scalar* PSIxR_elemwise = new scalar[Ne_AUG*(D+N_N)*M_G*N_s]; mem_counter.addToCPUCounter(Ne_AUG*(D+N_N)*M_G*N_s*sizeof(scalar));
  scalar* PSIxR_biased_global = new scalar[M_T*2*M_G*N_reco]; mem_counter.addToCPUCounter(M_T*2*M_G*N_reco*sizeof(scalar));
  scalar* PSIxR_biased_elemwise = new scalar[Ne_AUG*(D+N_N)*2*M_G*N_s]; mem_counter.addToCPUCounter(Ne_AUG*(D+N_N)*2*M_G*N_s*sizeof(scalar));
  //10/23/2017: ADding specific face-average operation for HAG scheme
  scalar* FaceAvg_global = new scalar[M_T*M_G*N_reco]; mem_counter.addToCPUCounter(M_T*M_G*N_reco*sizeof(scalar));
  for (int j = 0; j < Ne_AUG*(D+N_N)*M_G*N_s; j++){
    PSIxR_elemwise[j] = 0.0; }
  for (int j = 0; j < 2*Ne_AUG*(D+N_N)*M_G*N_s; j++){
    PSIxR_biased_elemwise[j] = 0.0; }


  //ideally, RecoPair would only need N_N entries per element, because it tells an
  //element who its interface in [0,M_T) is for a given local face index in [0,N_N).
  //Unfortunately, some of the global interfaces are duplicates of eachother,
  //so I need more than N_N entries per element. The worst possible case is D duplicates
  //touching one element, so I'm using N_N+D. PSIxR_elemwise must also adhere to this system
  int* RecoPair = new int[Ne_AUG*(D+N_N)]; mem_counter.addToCPUCounter(Ne_AUG*(D+N_N)*sizeof(int));
  for (int e = 0; e < Ne_AUG; e++)
    {
      //for (int s = 0; s < 2+N_N; s++)
      for (int s = 0; s < D + N_N; s++)
	{
	  RecoPair[e*(D+N_N) + s] = -1;
	}
    }
 
  
  int sequential_FaceFromElem[Ne_AUG];// = new int[N_E];
  for (int e = 0; e < Ne_AUG; e++){
    //    for (int s = 0; s < N_N; s++){
    sequential_FaceFromElem[e] = -1; }
  

//No need for a barrier; the only communication necessary
//for recovery is the XYZNodes_extended matrix,
//which is guaranteed to already be populated because
//of the Waitall command at end of CommunicateGhosts routine

   int sIndex_forReco[M_T*2];// = new int[M_T*2];
   for (int t = 0; t < M_T; t++)
     {
       sIndex_forReco[t*2 + 0] = 0;
       sIndex_forReco[t*2 + 1] = 0;
     }
   
	    //return 0;

   int PerformRecovery = 1;
   
   //PerformRecovery = 0; //10/10/2017: Setting this to zero because recovery is not ready for boundary conditions
   if (PerformRecovery == 1)
     {

       //Build some common entities to be used by all recovery routines.
       //The 2Sh stands for "to shell"
       fullMatrix<double> points_2Sh, weight_2Sh; //quadrature info on reference element
       if     (inputs.getElemType() == "lin") gaussIntegration::getLine(20, points_2Sh, weight_2Sh);
       else if(inputs.getElemType() == "tri") gaussIntegration::getTriangle(15, points_2Sh, weight_2Sh);
       else if(inputs.getElemType() == "qua") gaussIntegration::getQuad(15, points_2Sh, weight_2Sh);
       else if(inputs.getElemType() == "hex") gaussIntegration::getHexahedron(4*order+1, points_2Sh, weight_2Sh);
       int N_superG = points_2Sh.size1();
       fullMatrix<scalar> phiRef_2Sh(N_superG, N_s);
       fullMatrix<double> phiD_2Sh(N_superG, N_s);
       fullMatrix<scalar> dphi_2Sh(N_superG*D , N_s);
       basis->f (points_2Sh, phiD_2Sh);
       for(int g = 0; g < N_superG; g++){
	 for(int i = 0; i < N_s; i++){
	   phiRef_2Sh(g,i) = (scalar)phiD_2Sh(g,i);
	 }
       }   
  
       double grads_2Sh[N_s][3];  
       for(int g = 0; g < N_superG; g++){
	 basis->df(points_2Sh(g,0),points_2Sh(g,1),points_2Sh(g,2),grads_2Sh);
	 for(int alpha = 0; alpha < D; alpha ++){
	   for(int i = 0; i < N_s; i++){
	     dphi_2Sh(g*D+alpha,i) = (scalar)grads_2Sh[i][alpha];  //see paper for indexing p.6
	   }	  
	 }    
       }



       printf("pr=%d: Beginning PSIxR treatment for the full set of interfaces, M_T=%d\n",myid,M_T);
       //10/10/2017: Need to finally adjust this procedure for handling boundary interfaces
       for (int t = 0; t < M_T; t++)
       // for (int t = 0; t < 1; t++)
	 {
	   if (verbose > 0) {printf("\n========PSIxR treatment, pr=%d, interface %d==========\n",myid, t);}
	   const simpleInterface &face = interfaces[t]; //fetches the interface corresponing to index t
	   int elIndex[2];
	   /*
	   //    printf("face = %d\n",t);
	   for (int side = 0; side < 2; side++)
	   {
	   const simpleElement *el = face.getElement(side); //grabs the element on a specific side of the interface
	   elIndex[side] = ElementMap.at(el->getId()); //gets the element's index
	   }
       
	   //Now: Determine where PSIxR_elemwise should place this interface's sub-structure
	   //within the local element's stuctutesre
	   int omA = elIndex[0];
	   int omB = elIndex[1];
	   */
	   int omA = BR2_Map[t*4 + 0];
	   int omB = BR2_Map[t*4 + 2];
	   elIndex[0] = omA;
	   elIndex[1] = omB;
	   int BoundaryTag = 0;
	   if (omA == omB)
	     {
	       BoundaryTag = 1;
	       if (verbose > 0)
		 {
		   printf("In PSIxR treatment: Discovered that interface %d is a boundary interface, sitting on element omA=%d\n", t, omA); 
		 }
	     }
	   //10/10/2017: of omA==omB, then we are on a boundary interface.
	   //Instead of standard recovery, set the recovery operator
	   //to get the element's trace along the boundary

	   //The normal out of the interface
	   scalar* FaceNormal = new scalar[D];
	   for (int a = 0; a < D; a++)
	     {
	       FaceNormal[a] = normals(a,t);
	     }
      
	   /*
	     int HA = -1;//tells the element which of its N_N faces match the global interface t
	     int HB = -1;
	     //for (int H = 0; H < N_N; H++)
	     for (int H = 0; H < N_N; H++)
	     {
	     if (Alt_FaceFromElem[omA*(N_N) + H] == t)
	     {
	     HA = H;
	     }
	     if (Alt_FaceFromElem[omB*(N_N) + H] == t)
	     {
	     HB = H;
	     }
	     }
	     if (HA == -1)
	     {
	     printf("\nCATASTROPHE!!!HA remains -1 in PSIxR treatments\n");
	     }
	     if (verbose == 1)
	     {
	     if (HB == -1)
	     {
	     printf("\nHB remains -1 in PSIxR treatment\n");
	     }
	     }
	     //In periodic case, some interfaces are double-numbered. So, if HA or HB is negative 1,
	     //it means one of the elements does not know that it borders the interface t.
	     int duplicate_yes = 0;
	     int duplicate_index = 2*M_T;
	     if (HB == -1)
	     {
	     duplicate_yes = 1;
	     int FOUND = 0;
	     if (verbose == 1){printf("Searching for the duplicate interface to properly set HB\n");}
	     for (int t2 = 0; t2 < M_T; t2++)
	     {
	     if (t != t2 && FOUND == 0)
	     {
	     if ((BR2_Map[t2*4 + 0] == omA && BR2_Map[t2*4 + 2] == omB) || (BR2_Map[t2*4 + 0] == omB && BR2_Map[t2*4 + 2] == omA))
	     {
	     //we found the duplicate interface! Now, find the appropriate HB value 
	     //corresponding to this interface
	     for (int H = 0; H < N_N; H++)
	     {
	     if (Alt_FaceFromElem[omB*N_N + H] == t2)
	     {
	     HB = H;
	     if (verbose == 1) {printf("Found duplicate interface at t2=%d, HB = %d\n", t2, HB);}
	     FOUND = 1;
	     duplicate_index = t2;
	     }
	     }
	     }
	     }
	     }

	     }
	   */
	

	   //Let each element known how many of its interfaces have been treated
   
	   sequential_FaceFromElem[elIndex[0]] += 1;
	   sequential_FaceFromElem[elIndex[1]] += 1;
	   //printf("pr=%d, sequential_FaceFromElem[%d (out of Ne_AUG=%d)] incremented to %d\n",myid, elIndex[0], Ne_AUG, sequential_FaceFromElem[elIndex[0]]);
	   //printf("pr=%d, sequential_FaceFromElem[%d (out of Ne_AUG=%d)] incremented to %d\n",myid, elIndex[1], Ne_AUG, sequential_FaceFromElem[elIndex[1]]);

	   int sIndex[2];
	   //sIndex[0] = HA;
	   //sIndex[1] = HB;

	   sIndex[0] = sequential_FaceFromElem[elIndex[0]];
	   sIndex[1] = sequential_FaceFromElem[elIndex[1]];

	   RecoPair[elIndex[0]*(D+N_N) + sIndex[0]] = t;
	   RecoPair[elIndex[1]*(D+N_N) + sIndex[1]] = t;
	   if (verbose == 1)
	     {
	       printf("RecoPair[%d] = %d\n", elIndex[0]*(D+N_N) + sIndex[0] , t);
	       printf("RecoPair[%d] = %d\n", elIndex[1]*(D+N_N) + sIndex[1] , t);
	       printf("RecoPair(%d,%d) = %d\n", elIndex[0] , sIndex[0] , t);
	       printf("RecoPair(%d,%d) = %d\n", elIndex[1] , sIndex[1] , t);
	       //    printf("t=%d, omA=%d, omB=%d: elA storing contribution is substructure %d, elB storing contribution in substructure %d\n", t, omA, omB, HA, HB);
	     }

	   if (sIndex[0] >= N_N+D || sIndex[1] >= N_N+D)
	     {
	       printf("CATASTROPHE in building recovery operators\n");
	       printf("sIndex0=%d, sIndex1=%d, N_N+D=%d\n", sIndex[0], sIndex[1], N_N+D); fflush(stdout);
	       return 0;
	       exit(1);
	     }
	   


	   //sIndex is in [0,N_N); it tells the code where to put
	   //the partial recovery operator in the PSIxR_elemwise matrix

	   //      printf("Element A = %d, Element B = %d\n", elIndex[0], elIndex[1]);
	   //Necessities for calling GetPSUxR_oneface:
	   //solution node geometry, element A
	   //solution node geometry, element B
	   //quadrature node locations, face
	   //     printf("Preparing to populate xGeo_A, xGeo_B, xGQ_Face\n");
	   scalar* xGeo_A = new scalar[N_s*D];
	   scalar* xGeo_B = new scalar[N_s*D];
	   scalar* xGQ_Face = new scalar[M_G*D];
	
	   //    printf("Allocated room for xGeo_A...xGQ_face\n");
	   for (int k = 0; k < N_s; k++)
	     {
	       for (int a = 0; a < D; a++)
		 {
		   //	      printf("XYZnodes(%d, %d) = %f\n", k, elIndex[0]*D + a, XYZNodes(k, elIndex[0]*D + a));
		   //    printf("XYZnodes(%d, %d) = %f\n", k, elIndex[1]*D + a, XYZNodes(k, elIndex[1]*D + a));
		   xGeo_A[k*D + a] = XYZNodes(k, elIndex[0]*D + a);
		   xGeo_B[k*D + a] = XYZNodes(k, elIndex[1]*D + a);
		 }
	     }
	   //PEJ 11/28/2017: if it's a ghost element (off-partition), populate
	   //xGeo_B based on the corresponding flesh element, whose geometry
	   //is stored in XYZNodes_extended_Matrix
	   if (omB >= N_E)
	     {
	       for (int k = 0; k < N_s; k++){
		 for (int a = 0; a < D; a++){
		   xGeo_B[k*D + a] = XYZNodes_extended_Matrix(k, elIndex[1]*D + a);
		   //xGeo_B[k*D + a] = XYZNodes_extended[elIndex[1]*D*N_s + a*N_s + k]; 
		 }
	       }
	     }
	   //To adjust for periodicity: First must identify if the edge is periodic
	   //    printf("Finished populating xGeo_A and xGeo_B\n");
	   for (int g = 0; g < M_G; g++)
	     {
	       for (int a = 0; a < D; a++)
		 {
		   //I think the first half of columns per interface of XYGF is from element 0, second column is quadrature nodes in reverse direction
		   xGQ_Face[g*D + a] = XYZGF(g, t*2*D + a);
		 }
	     }
	   /*
	   //To adjust for periodicity: First must identify if the edge is periodic
	   scalar minFaceToB[D];
	   for (int a = 0; a < D; a++)
	   {
	   minFaceToB[a]= 9001;
	   }
	   //in each coordinate direction, get distance between first face quadratude node and centroid of element B
	   for (int k = 0; k < N_s; k++)
	   {
	  
	   }
	   */
	   /*
	     printf("Solution nodes, element A:\n");
	     for (int k = 0; k < N_s; k++)
	     {
	     printf("Node %d: (x,y,z) = (",k);
	     for (int a = 0; a < D; a++)
	     {
	     printf("%f, ", xGeo_A[k*D + a]);
	     }
	     printf(")\n");
	     }
	     printf("Solution nodes, element B:\n");
	     for (int k = 0; k < N_s; k++)
	     {
	     printf("Node %d: (x,y,z) = (",k);
	     for (int a = 0; a < D; a++)
	     {
	     printf("%f, ", xGeo_B[k*D + a]);
	     }
	     printf(")\n");
	     }
	     printf("Quadrature nodes, interface:\n");
	     for (int g = 0; g < M_G; g++)
	     {
	     printf("Node %d: (x,y,z) = (",g);
	     for (int a = 0; a < D; a++)
	     {
	     printf("%f, ", xGQ_Face[g*D + a]);
	     }
	     printf(")\n");
	     }
	   */
	   //storage scheme:
	   //PSIxR_elemwise(element index, local face index, face quadrature point, contributing DOF)


	   timers.start_timer(54);
	   scalar* PSIxR_local = new scalar[M_G*N_reco];
	   scalar* PSIxR_A_local = new scalar[M_G*N_s];
	   scalar* PSIxR_B_local = new scalar[M_G*N_s];

	   //PEJ 11/7/2017: Detailed deJ arrays of the two 
	   
	   PSIxR_shell(inputs, timers, N_N, N_s, order,/* m,*/ M_G, GradMethod_Vis, N_superG, phiRef_2Sh, dphi_2Sh, points_2Sh, weight_2Sh, xGeo_A, xGeo_B, xGQ_Face, FaceNormal, BoundaryTag, psi, Better_InvMap, t, M_s, omA, omB, PSIxR_local, PSIxR_A_local, PSIxR_B_local);
	   //GetPSIxR_oneface(inputs, N_s, order, N_reco, M_G, GradMethod_Vis, xGeo_A, xGeo_B, xGQ_Face, FaceNormal, BoundaryTag, psi, Better_InvMap, t, M_s, omA, omB, PSIxR_local, PSIxR_A_local, PSIxR_B_local);
	   //GetPSIxR_oneface(inputs, N_s, order, N_reco, M_G, xGeo_A, xGeo_B, xGQ_Face, FaceNormal, PSIxR_local, PSIxR_A_local, PSIxR_B_local);
	   for (int row = 0; row < M_G; row++)
	     {
	       for (int col = 0; col < 2*N_s; col++)
		 {
		   PSIxR_global[t*(M_G*2*N_s) + row*(2*N_s) + col] = PSIxR_local[row*2*N_s + col];
		 }
	       for (int col = 0; col < N_s; col++)
		 {
		   PSIxR_elemwise[elIndex[0]*((D+N_N)*M_G*N_s) + sIndex[0]*(M_G*N_s) + row*N_s + col] = PSIxR_local[row*2*N_s + col];
		   PSIxR_elemwise[elIndex[1]*((D+N_N)*M_G*N_s) + sIndex[1]*(M_G*N_s) + row*N_s + col] = PSIxR_local[row*2*N_s + N_s + col];
		   /*
		     PSIxR_elemwise[elIndex[0]*(N_N*M_G*N_s) + sIndex[0]*(M_G*N_s) + row*N_s + col] = PSIxR_A_local[row*N_s + col];
		     PSIxR_elemwise[elIndex[1]*(N_N*M_G*N_s) + sIndex[1]*(M_G*N_s) + row*N_s + col] = PSIxR_B_local[row*N_s + col];
		   */
		 }
	     }
   

      
	   //Repeat last few lines again to get a strict global FaceAvg operator; this is for
	   //HAG scheme, where I want interface gradient to use average but element
	   //gradient to use recovery; 1 is strict average in PSIxR_shell
	   if (HAG == 1)
	     {
	       PSIxR_shell(inputs, timers, N_N, N_s, order, /*m,*/ M_G, 1, N_superG, phiRef_2Sh, dphi_2Sh, points_2Sh, weight_2Sh, xGeo_A, xGeo_B, xGQ_Face, FaceNormal, BoundaryTag, psi, Better_InvMap, t, M_s, omA, omB, PSIxR_local, PSIxR_A_local, PSIxR_B_local);
	       for (int row = 0; row < M_G; row++)
		 {
		   for (int col = 0; col < 2*N_s; col++)
		     {
		       FaceAvg_global[t*(M_G*2*N_s) + row*(2*N_s) + col] = PSIxR_local[row*2*N_s + col];
		     }
		   //Don't need the elementwise matrix, just the global operator
		 }
	     }
	   delete[] PSIxR_local;
	   delete[] PSIxR_A_local;
	   delete[] PSIxR_B_local;
	   timers.stop_timer(54);

	   timers.start_timer(55);
	   //Now, the biased recovery operators
	   int N_EP;  //nodal equivalence points per element (may need to change for 3D)
	   if     (inputs.getElemType() == "lin") {N_EP=1;}
	   else if(inputs.getElemType() == "tri") {N_EP=order+1;}
	   else if(inputs.getElemType() == "qua") {N_EP=order+1;}
	   else if(inputs.getElemType() == "hex") {N_EP=pow(order+1,2);}
	   int N_icb = N_EP + N_s;
	   //N_icb = N_icb - 1;
	   if (verbose == 1) {printf("order=%d, N_icb=%d, bout to assemble biased recovery operator for t=%d\n",order, N_icb, t);}
	   scalar* PSIxR_biA_local = new scalar[M_G*2*N_s];
	   scalar* PSIxR_biA_A_local = new scalar[M_G*N_s];
	   scalar* PSIxR_biA_B_local = new scalar[M_G*N_s];
	   scalar* PSIxR_biB_local = new scalar[M_G*2*N_s];
	   scalar* PSIxR_biB_A_local = new scalar[M_G*N_s];
	   scalar* PSIxR_biB_B_local = new scalar[M_G*N_s];
	   //return 0;
	   //zero the PSIxR arrays and set nonzero only if using ICBN
	   for (int g = 0; g < M_G; g++){
	     for (int col = 0; col < N_s; col++){
	       PSIxR_biA_local[g*2*N_s + 0*N_s + col] = 0.0;
	       PSIxR_biA_local[g*2*N_s + 1*N_s + col] = 0.0;
	       PSIxR_biA_A_local[g*N_s + col] = 0.0;
	       PSIxR_biA_B_local[g*N_s + col] = 0.0;
	       PSIxR_biB_local[g*2*N_s + 0*N_s + col] = 0.0;
	       PSIxR_biB_local[g*2*N_s + 1*N_s + col] = 0.0;
	       PSIxR_biB_A_local[g*N_s + col] = 0.0;
	       PSIxR_biB_B_local[g*N_s + col] = 0.0; }}

#ifdef ICBN
	   //GetPSIxR_Biased_oneface_Nov2017(inputs, N_s, order, /*m,*/ N_icb, M_G, xGeo_A, xGeo_B, xGQ_Face, FaceNormal, PSIxR_biA_local, PSIxR_biB_local, PSIxR_biA_A_local, PSIxR_biA_B_local, PSIxR_biB_A_local, PSIxR_biB_B_local);
	   GetPSIxR_Biased_oneface_Nov2017(inputs, N_N, timers, N_s, order, /*m,*/ N_icb, M_G, basis, N_superG, phiRef_2Sh, dphi_2Sh, points_2Sh, weight_2Sh, xGeo_A, xGeo_B, xGQ_Face, FaceNormal, PSIxR_biA_local, PSIxR_biB_local, PSIxR_biA_A_local, PSIxR_biA_B_local, PSIxR_biB_A_local, PSIxR_biB_B_local);

#endif
	   //Put the biased recovery operator in to global storage.
	   //This one is twice the size as standard recovery operator
	   //because each element is part of 2 recovered solutions per interface.
	   //Storage scheme for those interested:
	   //PSIxR_biased_elemwise(element index, local face index, A/B dominance, face quadrature point, contributing DOF)
	   for (int row = 0; row < M_G; row++)
	     {
	       for (int dom = 0; dom < 2; dom++)
		 {
		   for (int col = 0; col < 2*N_s; col++)
		     {
		       PSIxR_biased_global[t*(2*M_G*2*N_s) + 0*M_G*2*N_s + row*(2*N_s) + col] = PSIxR_biA_local[row*2*N_s + col];
		       PSIxR_biased_global[t*(2*M_G*2*N_s) + 1*M_G*2*N_s + row*(2*N_s) + col] = PSIxR_biB_local[row*2*N_s + col];
		     }
		 }
	       for (int col = 0; col < N_s; col++)
		 {
		   //A-dominant component, element A storage:
		   PSIxR_biased_elemwise[elIndex[0]*((D+N_N)*2*M_G*N_s) + sIndex[0]*(2*M_G*N_s) + 0*(M_G*N_s) + row*N_s + col] = PSIxR_biA_A_local[row*N_s + col];
		   //B-dominant component, element A storage:
		   PSIxR_biased_elemwise[elIndex[0]*((D+N_N)*2*M_G*N_s) + sIndex[0]*(2*M_G*N_s) + 1*(M_G*N_s) + row*N_s + col] = PSIxR_biB_A_local[row*N_s + col];
		   //A-dominant component, element B storage:
		   PSIxR_biased_elemwise[elIndex[1]*((D+N_N)*2*M_G*N_s) + sIndex[1]*(2*M_G*N_s) + 0*(M_G*N_s) + row*N_s + col] = PSIxR_biA_B_local[row*N_s + col];
		   //B-dominant component, element B storage:
		   PSIxR_biased_elemwise[elIndex[1]*((D+N_N)*2*M_G*N_s) + sIndex[1]*(2*M_G*N_s) + 1*(M_G*N_s) + row*N_s + col] = PSIxR_biB_B_local[row*N_s + col];
		 }
	     }
	   delete[] PSIxR_biA_local;
	   delete[] PSIxR_biA_A_local;
	   delete[] PSIxR_biA_B_local;
	   delete[] PSIxR_biB_local;
	   delete[] PSIxR_biB_A_local;
	   delete[] PSIxR_biB_B_local;

	   timers.stop_timer(55);
	   
	   delete[] FaceNormal;
	   delete[] xGeo_A;
	   delete[] xGeo_B;
	   delete[] xGQ_Face;
	   //return 0;

	 }
       printf("pr=%d: Passed the interface loop for PSIxR and PSIxR_biased\n",myid);
       //return 0;
     }

if (verbose == 1) 
  {
    printf("Global PSIxR:\n");
    for (int t = 0; t < M_T; t++)
      // for (int t = 0; t < 1; t++)
      {
	int omA = BR2_Map[t*4+0];
	int omB = BR2_Map[t*4+2];
	printf("pr=%d, Interface %d: partition's element indices/N_E are %d, %d:\n",myid, t, omA, omB);
	for (int row = 0; row < M_G; row++)
	  {
	    printf("Quadrature node %d: ",row);
	    for (int col = 0; col < 2*N_s; col++)
	      {
		printf("%f, ",PSIxR_global[t*M_G*2*N_s + row*2*N_s + col]);
	      }
	    printf("\n");
	  }
      }
  }

//return 0;
     
     
     
     /*
       Recovery matrices are looking good. Need to test.
     */
     
     if (verbose == 1 || RecoTest == 1)
       {
      printf("RecoPair before -1->0 fix. For use in initialization. :\n");
      for (int e = 0; e < Ne_AUG; e++)
	{
	  for (int s = 0; s < (D+N_N); s++)
	    {
	      printf("RecoPair(%d,%d) = %d\n", e, s, RecoPair[e*(D+N_N)+s]);
	    }
	}
    }
  //return 0;

  //Test the recovery procedure. This can be run on
  //all elements or just one. If the recovery procedure
  //gives erroneous values, it is likely the result
  //of an ill-conditioned matrix.
  scalar* omA_check = new scalar[N_s* N_F];
  scalar* omB_check = new scalar[N_s* N_F];
  scalar* RecoDOF = new scalar[N_F*2*N_s];
  scalar* UhCommon_test = new scalar[M_T*M_G*N_F];
  scalar* UicbA_test = new scalar[M_T*M_G*N_F];
  scalar* UicbB_test = new scalar[M_T*M_G*N_F];
  if (RecoTest == 1)
    {
      for (int row = 0; row < N_s; row++)
	{
	  for (int col = 0; col < N_F; col++)
	    {
	      omA_check[col*N_s + row] = 1.0;
	      omB_check[col*N_s + row] = 3.0;
	    }
	}
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  //U is in column-major form
	  RecoDOF[fc*2*N_s + k] = omA_check[fc*N_s + k]; }}
  
      //Element B DOF:
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  RecoDOF[fc*2*N_s + N_s + k] = omB_check[fc*N_s + k]; }}
      //Matrix-multiply for the recovered solution:
      for (int t = 0; t < M_T; t++){
	printf("\ninterface %d:\n",t);
	for (int g = 0; g < M_G; g++){ 
	  int row = t*M_G + g;
	  for (int fc = 0; fc < N_F; fc++){
	    //Initialize sum to zero:
	    scalar sum = 0.0;
	  
	    for (int k = 0; k < 2*N_s; k++){
	      //here, I am summing over the 2K DOF contained in the element union for the specific field
	      sum += PSIxR_global[t*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
	    } //end the k summation loop
	    //Relay the sum to global UhCommon storage:
	    UhCommon_test[row*N_F + fc] = sum;
	    printf("UhCommon_test(t=%d, g=%d, f=%d) = %f\n", t, g, fc, UhCommon_test[row*N_F + fc]);
	  }
	}
      }
      //Matrix-multiply for the biased recovered solution:
      for (int t = 0; t < M_T; t++){
	const simpleInterface &face = interfaces[t]; //fetches the interface corresponing to index t
	int elIndex[2];
	//    printf("face = %d\n",t);
	for (int side = 0; side < 2; side++)
	  {
	    const simpleElement *el = face.getElement(side); //grabs the element on a specific side of the interface
	    elIndex[side] = ElementMap.at(el->getId()); //gets the element's index
	  }
       
	//Now: Determine where PSIxR_elemwise should place this interface's sub-structure
	//within the local element's stuctutesre
	int omA = elIndex[0];
	int omB = elIndex[1];
	scalar nx = normals(0, t);
	printf("\ninterface %d: omA=%d, omB=%d, nx=%f\n",t,omA,omB,nx);
	for (int g = 0; g < M_G; g++){ 
	  int row = t*M_G + g;
	  for (int fc = 0; fc < N_F; fc++){
	    //Initialize sum to zero:
	    scalar sum = 0.0;

	    for (int k = 0; k < 2*N_s; k++){
	      //here, I am summing over the 2K DOF contained in the element union for the specific field
	      sum += PSIxR_biased_global[t*(2*M_G*2*N_s) + 0*M_G*2*N_s + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
	    } //end the k summation loop
	    //Relay the sum to global UhCommon storage:
	    UicbA_test[row*N_F + fc] = sum;
	    printf("UicbA_test(t=%d, g=%d, f=%d) = %f\n", t, g, fc, UicbA_test[row*N_F + fc]);
	    //Repear for B-dominant solution
	    sum = 0.0;

	    for (int k = 0; k < 2*N_s; k++){
	      //here, I am summing over the 2K DOF contained in the element union for the specific field
	      sum += PSIxR_biased_global[t*(2*M_G*2*N_s) + 1*M_G*2*N_s + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
	    } //end the k summation loop
	    //Relay the sum to global UhCommon storage:
	    UicbB_test[row*N_F + fc] = sum;
	    printf("UicbB_test(t=%d, g=%d, f=%d) = %f\n", t, g, fc, UicbB_test[row*N_F + fc]);

	  } //end field variable loop
	} //end face quadr. point loop
      } //end interface loop
    } //end of RecoTest == 1 case
  //return 0;
  delete[] omA_check;
  delete[] omB_check;



  //No more use for UicbA_test, UicbB_test, UhCommon_test
  
  delete[] UhCommon_test;
  delete[] UicbA_test;
  delete[] UicbB_test;

  //Can also get rid of PSIxR_biased_global. May want to use it in future for HAGAR scheme, though
  

  //return 0;
  //look at map and invmap
  /*
  printf("Inspecting invmap\n");
  for (int e = 0; e < N_E; e++)
    {
      printf("Element %d, solidx:\n",e);
      for (int i = 0; i < M_s*N_N; i++)
	{
	  printf("[%d]->[%d],  ",i,h_invmap[((e*N_F+0)*M_s*N_N + i)*2+0]);
	  //	  printf("[%d]->[%d],  ",i,h_invmap[((e*N_F+1)*M_s*N_N+i)*2+0]);
	}
    }
  for (int e = 0; e < N_E; e++)
    {
      printf("Element %d, qidx:\n",e);
      for (int i = 0; i < M_s*N_N; i++)
	{
	  printf("[%d]->[%d],  ",i,h_invmap[((e*N_F+1)*M_s*N_N+i)*2+1]);
	}
    }
  */
  /*
  //Put Alt_FaceFromElem_serial in matrix form
  fullMatrix<scalar> Alt_FaceFromElem (N_E, N_N);
  for (int e = 0; e < N_E; e++)
    {
      for (int H = 0; H < N_N; H++)
	{
	  Alt_FaceFromElem(e,H) = Alt_FaceFromElem_serial[e*N_N + H];
	}
    }
  */
  //return 0;

  /*
  printf("Marc's node-node mapping\n");
  for (int t = 0; t < M_T; t++)
    {
      printf("Interface %d\n",t);
      for (int j = 0; j < M_s; j++)
	{
	  for (int f = 0; f < N_F; f++)
	    {
	      printf("\tnode=%d,field=%d\n",j,f);
	      int face_node;
	      for (int d = 0; d < 2; d++)
		{
		  face_node = t*N_F*2*M_s + f*2*M_s + d*M_s + j;//   ((t*N_F+fc)*2+d)*M_s+j;
		  printf("\t\tSide %d: face_node = %d, map[face_node] = %d\n",d,face_node, h_map[face_node]);
		}
	    }
	}
      printf("\n");
    }
  */
  //END BR2 edit

  
  //printf("Preparing for data transformation from fullMatrix to pointer form\n"); fflush(stdout);
  //
  //  We need to transform the data in fullMatrix to a pointer form to
  //  transfer to GPU note: it's got to be column major sorted
  //
  scalar* h_phi     = new scalar[N_G*N_s];       makeZero(h_phi,N_G*N_s);       mem_counter.addToCPUCounter(N_G*N_s*sizeof(scalar));
  scalar* h_phi_w   = new scalar[N_G*N_s];       makeZero(h_phi_w,N_G*N_s);     mem_counter.addToCPUCounter(N_G*N_s*sizeof(scalar));    
  scalar* h_dphi    = new scalar[D*N_G*N_s];     makeZero(h_dphi,D*N_G*N_s);	mem_counter.addToCPUCounter(D*N_G*N_s*sizeof(scalar)); 
  scalar* h_dphi_w  = new scalar[D*N_G*N_s];     makeZero(h_dphi_w,D*N_G*N_s);  mem_counter.addToCPUCounter(D*N_G*N_s*sizeof(scalar));
//printf("Succeded with h_phi, h_phi_w, h_dphi, and h_dphi_w\n"); fflush(stdout);
  scalar* h_normals = new scalar[D*M_T];         makeZero(h_normals,D*M_T);     mem_counter.addToCPUCounter(D*M_T*sizeof(scalar)); 
  scalar* h_psi     = new scalar[M_G*M_s];       makeZero(h_psi,M_G*M_s);	mem_counter.addToCPUCounter(M_G*M_s*sizeof(scalar)); 
//printf("Check 1\n"); fflush(stdout);
  scalar* h_psi_w   = new scalar[M_G*M_s];       makeZero(h_psi_w,M_G*M_s);	mem_counter.addToCPUCounter(M_G*M_s*sizeof(scalar)); 
//printf("Check 2\n"); fflush(stdout);
  scalar* h_J       = new scalar[N_E];           makeZero(h_J,N_E);             mem_counter.addToCPUCounter(N_E*sizeof(scalar)); 
// printf("Check 3\n"); fflush(stdout);
  //scalar* h_invJac  = new scalar[N_G*D*N_E*D];   makeZero(h_invJac,N_G*D*N_E*D);mem_counter.addToCPUCounter(N_G*D*N_E*D*sizeof(scalar)); 
  scalar* h_invJac  = new scalar[N_G*D*Ne_AUG*D];   makeZero(h_invJac,N_G*D*Ne_AUG*D);mem_counter.addToCPUCounter(N_G*D*Ne_AUG*D*sizeof(scalar)); 
// printf("Check 4\n"); fflush(stdout);
  scalar* h_JF      = new scalar[2*M_T];         makeZero(h_JF, 2*M_T);         mem_counter.addToCPUCounter(2*M_T*sizeof(scalar)); 
// printf("Check 5\n"); fflush(stdout);
  
printf("pr=%d: Succeded with h_psi, h_psi_w, h_J, h_invJac, h_JF, and h_normals\n",myid);
#ifdef USE_CPU
  // Normal allocation of U
  scalar* h_U       = new scalar[N_s*(N_E+N_ghosts)*N_F];   makeZero(h_U,N_s*(N_E+N_ghosts)*N_F); mem_counter.addToCPUCounter(N_s*(N_E+N_ghosts)*N_F*sizeof(scalar));
#elif USE_GPU
  // Pinned allocation of U. Faster host <-> device transfer
  // http://devblogs.nvidia.com/parallelforall/how-optimize-data-transfers-cuda-cc/
  // If I run into an issue: cudaHostAlloc and cudaMallocHost ARE NOT THE SAME (http://stackoverflow.com/questions/8591577/does-cudafreehost-care-what-device-is-active-when-cudamallochost-is-used-to-allo)
  scalar* h_U;
  checkCuda(cudaMallocHost((void**)&h_U, N_s*(N_E+N_ghosts)*N_F*sizeof(scalar)));
#endif

  timers.start_timer(58);
  //No hang

  //return 0;
  // copy from the fullMatrix to the pointer format (column major)
  phi.copyMatrixToPointer(h_phi);
  phi_w.copyMatrixToPointer(h_phi_w);
  dphi.copyMatrixToPointer(h_dphi);
  dphi_w.copyMatrixToPointer(h_dphi_w);
  psi.copyMatrixToPointer(h_psi);
  psi_w.copyMatrixToPointer(h_psi_w);
  
  /*
  printf("h_psi organization:\n");
  for (int j = 0; j < M_G*M_s; j++)
    {
      printf("j=%d, psi=%f\n", j, h_psi[j]);
    }
  */
  //return 0;
  for(int e = 0; e < N_E; e++){
    h_J[e] = J(e,0);
  }

  //need h_invJac in ghost elements for mixed form appraoch
  //  for(int e = 0; e < N_E; e++){
  for(int e = 0; e < Ne_AUG; e++){
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
  //column-major-example
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      for(int i = 0; i < N_s; i++){
  	h_U[(e*N_F+fc)*N_s+i]= U(i,e*N_F+fc);
      }
    }
  }

  //PEJ 12/16/2017: Communicate the U vector
  //across partitions so I can 
  //test some operators
  if (RecoTest + MapTest + AuxSolveTest + SigFaceTest > 0)
    {
      communicator.CommunicateGhosts(N_F, h_U); //populate ghost elements from flesh elements
    }
  timers.stop_timer(58);

  /*
  printf("\nAlt_FaceFromElem:\n");
  for (int e = 0; e < N_E; e++)
    {
      printf("Element %d tGlo values are \n(", e);
      for (int H = 0; H < N_N; H++)
	{
	  printf("%d, ",Alt_FaceFromElem[e*N_N + H]);
	}
      printf(")\n");
    }
  */
  
  timers.start_timer(57);
  //PEJ edit 05/21/2017: Populate the mixed-form surface and volume contributions
  //for the auxiliary solve.
  //Need to copy Marc's style for declaring these large entities
  scalar* serial_MSigSurf = new scalar[N_E * (D*N_G) * (N_N*M_G)]; makeZero(serial_MSigSurf,N_E*D*N_G*N_N*M_G);     mem_counter.addToCPUCounter(N_E*D*N_G*N_N*M_G*sizeof(scalar)); 
  scalar* serial_MSigVol = new scalar[N_E * (D*N_G) * (N_s)];      makeZero(serial_MSigVol,N_E*D*N_G*N_s); mem_counter.addToCPUCounter(N_E*D*N_G*N_s*sizeof(scalar));
// printf("\npr=%d: Calling BuildSigMatrices\n",myid);
  //return 0;
  //BuildSigMatrices(inputs, N_E, order, N_G, M_G, N_N, M_s, h_invmap, pointsF, weightF, JF, XYZNodes, h_normals, h_invJac, J, FaceFromElem, BR2_Map, phi, h_Minv, serial_SigSurf, serial_SigVol);
  /*
  printf("JFace:\n");
  for (int t = 0; t < M_T; t++)
    {
      printf("tGlo=%d: JFace=%f\n",t, JF(t*2,0));
    }
  */
  BuildSigMatrices(inputs, QuadRule, N_E, N_s, order, N_G, M_G, N_N, M_s, h_invmap, pointsF, weightF, JF, XYZNodes, normals, h_invJac, J, Alt_FaceFromElem, BR2_Map, phi, h_Minv, psi_w, detJ_full, serial_MSigSurf, serial_MSigVol);
  printf("pr=%d: Finished with BuildSigMatrices\n",myid); fflush(stdout);

  //No Hang

  //return 0;
  if (verbose == 1)
    {
      printf("Global PSIxR:\n");
      for (int t = 0; t < 1/*M_T*/; t++)
	{
	  printf("Interface %d:\n",t);
	  for (int row = 0; row < M_G; row++)
	    {
	      printf("Quadrature node %d: ",row);
	      for (int col = 0; col < 2*N_s; col++)
		{
		  printf("%f, ",PSIxR_global[t*M_G*2*N_s + row*2*N_s + col]);
		}
	      printf("\n");
	    }
	}
    }

  //PEJ 10/10/2017: Account for Sigma calculation on boundary interfaces.
  /*
    Summary of strategy: Assuming I have Uhat for an element
    and UhCommon along its perimeter, I can use the summation combination of
    SurfForSigma and MColForSigma to get sigma coefficinets for the element.

    Then, use the face basis functions (psi in Marc's language, not recovery basis)
    to get the trace of the auxiliary variable along the boundary interface.
   */
  //serial_AuxHatSurf gets the auxiliary coefficient contributions from UhCommon distribution.
  //serial_AuxHatVol gets the auxiliary coefficient contributions from Uhat.
  scalar* serial_AuxHatSurf = new scalar[N_E * (D*N_s) * (N_N*M_G)]; //DELETE
  scalar* serial_AuxHatVol = new scalar[N_E * (D*N_s) * N_s]; //DELETE
  //With these populated, it's possible to get boundary gradient by multilyeing with shape functions.
  //This subroutine is like BuildSigMatrices except the multiplication by phi is ommitted
  BuildAuxHatMatrices(inputs, QuadRule, N_E, N_s, order, N_G, M_G, N_N, M_s, h_invmap, pointsF, weightF, JF, XYZNodes, normals, h_invJac, J, Alt_FaceFromElem, BR2_Map, phi, h_Minv, psi_w, detJ_full, serial_AuxHatSurf, serial_AuxHatVol);
  
  timers.stop_timer(57);

  fullMatrix<scalar> Get_Trace_From_Hat (M_G, N_s);
  scalar* serial_Uhat2GradBC = new scalar[M_B*D*M_G*N_s]; //DELETE
  scalar* serial_UhCommon2GradBC = new scalar[M_B*D*M_G*(N_N*M_G)]; //DELETE


  //No hang


  //Now, we must find a way to grab the boundary distribution of sigma
  for (int tB = 0; tB < M_B; tB++)
    {
      int tG = m.getBoundaryMap()[tB]; //global interface address
      //printf("tB=%d,tG=%d\n",tB,tG);
      int omA = BR2_Map[tG*4+0]; //supporting element
      if (verbose == 1) {printf("---pr=%d, tB=%d, tG=%d, omA=%d---\n",myid,tB,tG,omA);}
      for (int g = 0; g < M_G; g++){
	for (int k = 0; k < N_s; k++){
	  Get_Trace_From_Hat(g,k) = 0.0; }}
      
      for (int m = 0; m < M_s; m++) //supported node on boundary face
	{
	  int fc = 0;
	  int side = 0; //omA gets side 0
	  int face = ((tG*N_F+fc)*2+side)*M_s + m;
	  int Masterf0Index = Better_InvMap[face];
	  int row = Masterf0Index - omA*(N_F*N_s) - fc*N_s; //gives k = element index corresponding to the supported node
	  for (int g = 0; g < M_G; g++)
	    {
	      Get_Trace_From_Hat(g, row) = psi_w(g,m);
	      if (verbose == 1) {printf("Get_Trace_From_Hat(g=%d,k=%d)=%f\n",g,row,psi_w(g,m));}
	    }
	}
      //Now, boundaryGradient_from_Sighat is multiplied by the AuxHat matrices
      //to develop a structure that can directly calculate boundary gradient.
      for (int a = 0; a < D; a++) //gradient component
	{
	  for (int g = 0; g < M_G; g++) //quadrature point on the boundary interface
	    {
	      for (int k = 0; k < N_s; k++)
		{
		  scalar sum = 0;
		  for (int j = 0; j < N_s; j++)
		    {
		      sum += Get_Trace_From_Hat(g,j) * serial_AuxHatVol[omA*(D*N_s*N_s) + a*N_s*N_s + j*N_s + k];
		    }
		  serial_Uhat2GradBC[tB*D*M_G*N_s + a*M_G*N_s + g*N_s + k] = sum;
		  if (verbose == 1) {printf("serial_Uhat2GradBC[%d*D*M_G*N_s + %d*M_G*N_s + %d*N_s + %d]=%f\n",tB,a,g,k,sum);}
		}
	      for (int run = 0; run < N_N*M_G; run++)
		{
		  scalar sum = 0;
		  for (int j = 0; j < N_s; j++)
		    {
		      sum += Get_Trace_From_Hat(g,j) * serial_AuxHatSurf[omA*(D*N_s*N_N*M_G) + a*N_s*N_N*M_G + j*N_N*M_G + run];
		    }
		  serial_UhCommon2GradBC[tB*D*M_G*(N_N*M_G) + a*M_G*N_N*M_G + g*N_N*M_G + run] = sum;
		  if (verbose == 1) {printf("serial_UhCommon2GradBC[tB*D*M_G*(N_N*M_G) + a*M_G*N_N*M_G + g*N_N*M_G + run] = %f\n",sum);}
		}
	    }
	}
      //Check for maintainance of stationary state
      if (AuxSolveTest == 1)
	{
	  scalar FakeDOF[N_s];
	  scalar FakeUhCommon[N_N*M_G];
	  for (int k = 0; k < N_s; k++)
	    {
	      FakeDOF[k] = k*0.0;/*10.0;*/}
	  for (int j = 0; j < N_N*M_G; j++)
	    {
	      FakeUhCommon[j] = j*0.0;/*10.0;*/}
	  scalar FakeGrad[D*M_G];
	  for (int j = 0; j < D*M_G; j++)
	    {
	      FakeGrad[j] = 0.0;
	      scalar sum = 0;
	      for (int k = 0; k < N_s; k++)
		{
		  sum += serial_Uhat2GradBC[tB*D*M_G*N_s + j*N_s  + k] * FakeDOF[k];
		}
	      for (int run = 0; run < N_N*M_G; run++)
		{
		  sum += serial_UhCommon2GradBC[tB*D*M_G*(N_N*M_G)  + j*N_N*M_G + run] * FakeUhCommon[run];
		}
	      FakeGrad[j] = sum;
	      printf("Fake Grad(tb=%d, a*M_G+g=%d) = %f\n",tB,j,FakeGrad[j]);
	   
	    }
	}
    
      
    }
  
  //No Hang

  //  return 0;

  //PEJ 05/28/2017: Build interface sigma matrices
  //  printf("pr=%d: Declaring serial_SigFace_from_DOF\n",myid);
  timers.start_timer(56);
  scalar* serial_SigFace_from_DOF = new scalar[M_T*D*M_G*2*N_s]; makeZero(serial_SigFace_from_DOF,M_T*D*M_G*2*N_s); mem_counter.addToCPUCounter(M_T*D*M_G*2*N_s*sizeof(scalar));
  //For serial_sigFace: if not applying HAG scheme, pass in PSIxR_global
  if (HAG == 0)
    {
      SigFaceMatrices(M_T, M_G, N_s, N_G, N_N, M_s, Chi_Vis, phi, dphi, weight, J, psi, psi_w, normals, JF, Alt_FaceFromElem, BR2_Map, h_invmap, Better_InvMap, h_invJac, h_Minv, PSIxR_global, detJ_full, serial_SigFace_from_DOF);
    }
  else if (HAG == 1)
    {
      SigFaceMatrices(M_T, M_G, N_s, N_G, N_N, M_s, Chi_Vis, phi, dphi, weight, J, psi, psi_w, normals, JF, Alt_FaceFromElem, BR2_Map, h_invmap, Better_InvMap, h_invJac, h_Minv, FaceAvg_global, detJ_full, serial_SigFace_from_DOF);
    }
  else
    {
      printf("FAILURE: HAG out of bounds\n");
      exit(1);
    }

  timers.stop_timer(56);

  // printf("pr=%d: Ran SigFaceMatrices\n",myid);
  delete[] FaceAvg_global;
  //return 0;
  //No more need for PSIxR_Global, all that remains is to partition operators and run the simulation
  /*
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  //just for debug
#endif

  */
  //No Hang
  

  //07/02/2017 PEJ: splitting certain interface-based operations in to 
  //element-wise operations.
  int* BinarySideAddress = new int[Ne_AUG*(D+N_N)]; mem_counter.addToCPUCounter(Ne_AUG*(D+N_N)*sizeof(int));
  for (int j = 0; j < Ne_AUG*(D+N_N); j++) {
    BinarySideAddress[j] = 0; }
  //SigFace_from_Uhat_elemwise builds an element contrubution to send to the interface,
  //for each pair of elements sharing an interface. Thus, must include ghosts.
  scalar* SigFace_from_Uhat_elemwise = new scalar[Ne_AUG*(D+N_N)*D*M_G*N_s]; makeZero(SigFace_from_Uhat_elemwise, Ne_AUG*(D+N_N)*D*M_G*N_s); mem_counter.addToCPUCounter(Ne_AUG*(D+N_N)*D*M_G*N_s*sizeof(scalar));
  //zero the element-wise SigFace operation to be safe:
  for (int j = 0; j < Ne_AUG*(D+N_N)*D*M_G*N_s; j++) {
    SigFace_from_Uhat_elemwise[j] = 0.0; }
  //Re-zero the interface counter within each element
  for (int e = 0; e < Ne_AUG; e++){
    //    for (int s = 0; s < N_N; s++){
    sequential_FaceFromElem[e] = -1; }
  for (int t = 0; t < M_T; t++) //global interface address
    {
      const simpleInterface &face = interfaces[t]; //fetches the interface corresponing to index t
      int elIndex[2];
      //    printf("face = %d\n",t);
      
      //omB could be a ghost, my approach is build to handle it.
      int omA = BR2_Map[t*4 + 0];
      int omB = BR2_Map[t*4 + 2];
      elIndex[0] = omA;
      elIndex[1] = omB;

      if (verbose == 1){printf("pr=%d: Splitting SigFace_from_Uhat, interface %d, omA=%d, omB=%d\n",myid, t, omA, omB);}
      
      //increment each element's interface counter
      sequential_FaceFromElem[elIndex[0]] += 1;
      sequential_FaceFromElem[elIndex[1]] += 1;
       
       int sIndex[2];
       //sIndex[0] = HA;
       //sIndex[1] = HB;
       sIndex[0] = sequential_FaceFromElem[elIndex[0]];
       sIndex[1] = sequential_FaceFromElem[elIndex[1]];

       //BinarySideAddress tells an element for a given face, if it is elA or elB.
       BinarySideAddress[elIndex[0]*(D+N_N) + sIndex[0]] = 0;
       BinarySideAddress[elIndex[1]*(D+N_N) + sIndex[1]] = 1;

       //Send the SigFace operation to element-wise storage:
       //The indexing needs to match ordering used in dg/kernels_phil.cu/Uhat_to_GradCommon_v2.
       for (int a = 0; a < D; a++)
	 {
	   for (int g = 0; g < M_G; g++)
	     {
	       for (int k = 0; k < N_s; k++)
		 {
		   //the first element's contribution:
		   SigFace_from_Uhat_elemwise[elIndex[0]*(D+N_N)*M_G*D*N_s + 
					      sIndex[0] *(M_G*D*N_s) + 
					      g         *D*N_s + 
					      a         *N_s + 
					      k] = 
		     serial_SigFace_from_DOF[t*(D*M_G*2*N_s) + 
					     a*(M_G*2*N_s) + 
					     g*(2*N_s) +
					     k];
		   //The second element's contribution (extra N_s argument):
		   SigFace_from_Uhat_elemwise[elIndex[1]*(D+N_N)*M_G*D*N_s + 
					      sIndex[1] *(M_G*D*N_s) + 
					      g         *D*N_s + 
					      a         *N_s + 
					      k] = 
		     serial_SigFace_from_DOF[t*(D*M_G*2*N_s) + 
					     a*(M_G*2*N_s) + 
					     g*(2*N_s) + N_s + 
				     k];
		 }
	     }
	 } 
    }

 
  //Re-organize the interface-oriented SigFace operator for hopefully better looping
  scalar* serial_SigFace_from_Uhat_ROG = new scalar[M_T*D*M_G*2*N_s]; makeZero(serial_SigFace_from_Uhat_ROG,M_T*D*M_G*2*N_s); mem_counter.addToCPUCounter(M_T*D*M_G*2*N_s*sizeof(scalar));
  for (int t = 0; t < M_T; t++){
    for (int a = 0; a < D; a++){
      for (int g = 0; g < M_G; g++){
	for (int k = 0; k < 2*N_s; k++){
	  //still need to pass this thing in to DGsolver
	  serial_SigFace_from_Uhat_ROG[t*(M_G*D*2*N_s) + g*(D*2*N_s) + a*(2*N_s) + k] = serial_SigFace_from_DOF[t*(D*M_G*2*N_s) + a*(M_G*2*N_s) + g*(2*N_s) + k];
	}}}}

  //Do the biased recovery operation with element-specific storage
  scalar* Uicb_test = new scalar[2*M_T*M_G*N_F];
  scalar* UicbHalf_test = new scalar[2*2*M_T*M_G*N_F];
  scalar* Uhat_fake = new scalar[Ne_AUG*N_s*N_F];
  if (RecoTest == 1)
    {
      for (int e = 0; e < Ne_AUG; e++){
	for (int k = 0; k < N_s; k++){
	  for (int fc = 0; fc < N_F; fc++){
	    //Uhat_fake[(e*N_F+fc)*N_s + k] = e + 0.0; }}}
	    Uhat_fake[(e*N_F+fc)*N_s + k] = 2*((e+1) % 2) + 1.0; //either 1 or 3
	  }}}
      for (int j = 0; j < 2*2*M_T*M_G*N_F; j++){
	UicbHalf_test[j] = 0.0; }
      for (int j = 0; j < 2*M_T*M_G*N_F; j++){
	Uicb_test[j] = 0.0; }
      for (int e = 0; e < Ne_AUG; e++)
	{
	  scalar sum;
	  for (int s = 0; s < D+N_N; s++)
	    {
	      int tGlo = RecoPair[e*(D+N_N) + s];
	      int ArgSide = BinarySideAddress[e*(D+N_N) + s];
	      for (int g = 0; g < M_G; g++)
		{
		  for (int dom = 0; dom < 2; dom++)
		    {
		      for (int fc = 0; fc < N_F; fc++)
			{
			  sum = 0.0;
			  for (int k = 0; k < N_s; k++)
			    {
			      sum += PSIxR_biased_elemwise[e*((D+N_N)*2*M_G*N_s) + s*(2*M_G*N_s) + dom*M_G*N_s + g*N_s + k] * Uhat_fake[(e*N_F+fc)*N_s+k];
			    }
			  UicbHalf_test[tGlo*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + ArgSide] += sum;
			} //end fc loop
		    } //end dom loop
		} //end g loop
	    }
	}
      printf("Test of Uicb assembly:\n");
      for (int t = 0; t < M_T; t++) //global interface address
	{
	  const simpleInterface &face = interfaces[t]; //fetches the interface corresponing to index t
	  int elIndex[2];
	  //    printf("face = %d\n",t);
	  /*
	  for (int side = 0; side < 2; side++)
	    {
	      const simpleElement *el = face.getElement(side); //grabs the element on a specific side of the interface
	      elIndex[side] = ElementMap.at(el->getId()); //gets the element's index
	    }
       
	  //Now: Determine where PSIxR_elemwise should place this interface's sub-structure
	  //within the local element's stuctutesre
	  int omA = elIndex[0];
	  int omB = elIndex[1];
	  */
	   //omB could be a ghost, my approach is build to handle it.
	  int omA = BR2_Map[t*4 + 0];
	  int omB = BR2_Map[t*4 + 2];
	  elIndex[0] = omA;
	  elIndex[1] = omB;
	  scalar nx = normals(0, t);
	  printf("\npr=%d, interface %d: omA=%d, omB=%d, nx=%f\n",myid,t,omA,omB,nx);
	  for (int fc = 0; fc < N_F; fc++) //field
	    {
	      for (int dom = 0; dom < 2; dom++) //the side of interface whose ICB solution is being fully summed
		{
		  for (int g = 0; g < M_G; g++) //quadrature node on the interface
		    {
		      //zero and one are for the contributing elements
		      //Uicb_test[t*N_F*2*M_G + fc*2*M_G + dom*M_G + g] = UicbHalf_test[t*M_G*2*N_F*2 + g*2*N_F*2 + 0*N_F*2 + fc*2 + dom] + UicbHalf_test[t*M_G*2*N_F*2 + g*2*N_F*2 + 1*N_F*2 + fc*2 + dom];
		      Uicb_test[t*N_F*2*M_G + fc*2*M_G + dom*M_G + g] = UicbHalf_test[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 0] + UicbHalf_test[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 1];
		      printf("pr=%d,t=%d,f=%d,dom=%d,g=%d: left sum=%f, right sum=%f, Uicb=%f\n", myid,t, fc, dom, g, 
			     UicbHalf_test[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 0], 
			     UicbHalf_test[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 1], 
			     Uicb_test[t*N_F*2*M_G + fc*2*M_G + dom*M_G + g]);
		    } //end g loop
		} //end dom loop
	    } //end field variable loop
	} //end global interface loop

      //return 0;

      scalar* SigCommon_test = new scalar[M_T*D*M_G*N_F];
      //test the SigFace operator
      for (int row = 0; row < N_s; row++)
	{
	  for (int col = 0; col < N_F; col++)
	    {
	      omA_check[col*N_s + row] = 1.0;
	      omB_check[col*N_s + row] = 1.0;
	    }
	}
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  //U is in column-major form
	  RecoDOF[fc*2*N_s + k] = omA_check[fc*N_s + k]; }}
  
      //Element B DOF:
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  RecoDOF[fc*2*N_s + N_s + k] = omB_check[fc*N_s + k]; }}
  
      for (int t = 0; t < M_T; t++){
	printf("pr=%d, interface %d:\n",myid,t);
	for (int a = 0; a < D; a++) {
	  printf("\tfor a=%d:\n",a);
	  for (int g = 0; g < M_G; g++){ 
	    int row = t*D*M_G + a*M_G + g;
	    for (int fc = 0; fc < N_F; fc++){
	      //Initialize sum to zero:
	      scalar sum = 0.0;
	  
	      for (int k = 0; k < 2*N_s; k++){
		//here, I am summing over the 2K DOF contained in the element union for the specific field
		sum += serial_SigFace_from_DOF[t*(D*M_G*2*N_s) + a*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
	      } //end the k summation loop
	      //Relay the sum to global UhCommon storage:
	      SigCommon_test[row*N_F + fc] = sum;
	      printf("SigCommon_test(t=%d, a=%d, g=%d, f=%d) = %f\n", t, a, g, fc, SigCommon_test[row*N_F + fc]);
	    }
	  }
	}
      }
    } //end of RecoTest == 1 case
  
  
  
  //return 0;

  /*Here: for every (-1) in RecoPair, set it to zero.
    This means that in Uhat_to_UhCommon_v2 and Uhat_to_gradCommon_v2, the entries
    for tGlo=0 will he added to MANY times,
    but the additions will usually be zero
    because that is what PSIxR_elemwise (and SigFace_from_Uhat_elemwise)
    was initialized to, and RecoPair=-1 means that the PSIxR_elemwise entries are zero.*/
  for (int e = 0; e < Ne_AUG; e++)
    {
      for (int s = 0; s < (D+N_N); s++)
	{
	  if (RecoPair[e*(D+N_N) + s] == -1)
	    {
	      RecoPair[e*(D+N_N) + s] = 0;
	    }
	}
    }
  if (verbose == 1)
    {
      printf("pr=%d: RecoPair after -1->0 fix. For use in DG solver. :\n",myid);
      for (int e = 0; e < Ne_AUG; e++)
	{
	  for (int s = 0; s < (D+N_N); s++)
	    {
	      printf("pr=%d, RecoPair(%d,%d) = %d\n", myid, e, s, RecoPair[e*(D+N_N)+s]);
	    }
	}
    }


  
 
  //return 0;

  if (MapTest == 1)
    {
      printf("THE MAP TEST, pr=%d: Check if solution collocation to interface is effective\n",myid);
      //Test procedure for collocating the solution to the interfaces
      scalar* UF = new scalar[2*N_F*M_s*M_T];    makeZero(UF,2*N_F*M_s*M_T); 
      scalar* UintegF = new scalar[2*N_F*M_G*M_T];    makeZero(UintegF,2*N_F*M_G*M_T);
      // map U onto UF: requires Map, Ustar, UF and some integers for sizes, etc
      communicator.CommunicateGhosts(N_F, h_U); //populate ghost elements from flesh elements
      LmapToFace(M_s, M_T, N_s, h_map, h_U, UF); 
      blasGemm('N','N', M_G, M_T*N_F*2, M_s, 1, h_psi, M_G, UF, M_s, 0.0, UintegF, M_G);
      //LMap to face yields the nodal DOF of the face, NOT quadrature point solutions.
      //Then, the blasGemm call does matrix multiply to send solution to face quadrature points
      for (int t = 0; t < M_T; t++)
	{
	  printf("\t------pr=%d, Interface %d:------\n",myid,t);
	  const simpleInterface &face = interfaces[t]; //fetches the interface corresponing to index t
	  
	  int elIndex[2];
	  //    printf("face = %d\n",t);
	  /*
	  for (int side = 0; side < 2; side++)
	    {
	      const simpleElement *el = face.getElement(side); //grabs the element on a specific side of the interface
	      elIndex[side] = ElementMap.at(el->getId()); //gets the element's index
	    }
       
	  //Now: Determine where PSIxR_elemwise should place this interface's sub-structure
	  //within the local element's stuctutesre
	  int omA = elIndex[0];
	  int omB = elIndex[1];
	  */
	  int omA = BR2_Map[t*4+0];
	  int omB = BR2_Map[t*4+2];
	  elIndex[0] = omA;
	  elIndex[1] = omB;
	  int id1 = face.getClosureId(0);
	  int id2 = face.getClosureId(1);
	  const std::vector<int> &cl1 = closures[id1];
	  const std::vector<int> &cl2 = closures[id2];
#ifdef THREED
	  printf("\tpr=%d: The closure id for omA={%d,%d,%d,%d}, for omB={%d,%d,%d,%d}\n", myid, 
		 cl1[0], cl1[1], cl1[2], cl1[3],
		 cl2[0], cl2[1], cl2[2], cl2[3]);
#endif
	  for (int k = 0; k < M_s; k++)
	    {
	      //    printf("\t\t pr=%d: face node %d: omA_index=%d, omB_index=%d, rhoL=%f, rhoR=%f\n", myid, k, h_map[((t*N_F+0)*2+0)*M_s + k] - omA*N_F*N_s - 0*N_s, h_map[((t*N_F+0)*2+1)*M_s + k] - omB*N_F*N_s - 0*N_s, UF[t*N_F*2*M_s + 0*2*M_s + 0*M_s + k], UF[t*N_F*2*M_s + 0*2*M_s + 1*M_s + k]);
	    }
	  printf("\tpr=%d: The quadrature node solutions:\n", myid);
	  for (int g = 0; g < M_G; g++)
	    {
	      scalar diff = fabs(UintegF[t*N_F*2*M_G + 0*2*M_G + 0*M_G + g]-UintegF[t*N_F*2*M_G + 0*2*M_G + 1*M_G + g]);
	      if (diff < 0.0000001)
		{
		  //	  printf("\t\t pr=%d: quadrature node %d: rhoL=%f, rhoR=%f\n", myid, g, UintegF[t*N_F*2*M_G + 0*2*M_G + 0*M_G + g], UintegF[t*N_F*2*M_G + 0*2*M_G + 1*M_G + g]);
		}
	      else
		{
		  printf("\t\t pr=%d: quadrature node %d: rhoL=%f, rhoR=%f\t\t\t\t\tDISCREPANCY=%f\n", myid, g, UintegF[t*N_F*2*M_G + 0*2*M_G + 0*M_G + g], UintegF[t*N_F*2*M_G + 0*2*M_G + 1*M_G + g],diff);
		}
	    }
	}
      delete[] UintegF;
      delete[] UF;
    }
  
  // return 0;

  //With SigFace splitting and SigFace testing complete, we can delete serial_SigFace_from_DOF

  //testing phase concluded, delete.
  delete[] Uicb_test;
  delete[] UicbHalf_test;
  delete[] Uhat_fake;
  delete[] RecoDOF;

  //return 0;

  //Test the auxiliary variable solve:
  scalar* TestUhCommon = new scalar[M_T*M_G*N_F];
  scalar* TestAuxHat = new scalar[N_E*D*N_s*N_F];
  if (AuxSolveTest == 1)
    {
      printf("pr=%d, EXECUTING AUXILIARY SOLVE TEST\n",myid);
      for (int t = 0; t < M_T; t++)
	{
	  int omA = BR2_Map[t*4 + 2*0 + 0];
	  int omB = BR2_Map[t*4 + 2*1 + 0];
	  for (int g = 0; g < M_G; g++)
	    {
	      for (int fc = 0; fc < N_F; fc++)
		{
		  scalar sum = 0;
		  for (int j = 0; j < N_s; j++)
		    {
		      sum += PSIxR_global[t*(M_G*2*N_s) + g*(2*N_s) + j]       * h_U[(omA*N_F+fc)*N_s+j];
		      sum += PSIxR_global[t*(M_G*2*N_s) + g*(2*N_s) + N_s + j] * h_U[(omB*N_F+fc)*N_s+j];
		    }
		  TestUhCommon[t*M_G*N_F + g*N_F + fc] = sum;
		  //printf("Recovered Uh (t=%d, g=%d, fc=%d) (x,y)=(%f,%f) = %f\n",t,g,fc, XYZGF(g, t*2+0*D+0), XYZGF(g, t*2+0*D+1),TestUhCommon[t*M_G*N_F + g*N_F + fc]);
		}
	    }
	}
      //Correct boundary interfaces to 2.0
      for (int tB = 0; tB < M_B; tB++)
	{
	  //      int tG = m.boundaryMap()[tB];
	  int tG = m.getBoundaryMap()[tB]; //global interface address
	  for (int g = 0; g < M_G; g++)
	    {
	      for (int fc = 0; fc < N_F; fc++)
		{
		  TestUhCommon[tG*M_G*N_F + g*N_F + fc] = 2.0;
		}
	    }
	}
      for (int t = 0; t < M_T; t++)
	{
	  for (int g = 0; g < M_G; g++)
	    {
	      for (int fc = 0; fc < N_F; fc++)
		{
		  printf("Recovered Uh (t=%d, g=%d, fc=%d) (x,y)=(%f,%f) = %f\n",t,g,fc, XYZGF(g, t*2+0*D+0), XYZGF(g, t*2+0*D+1),TestUhCommon[t*M_G*N_F + g*N_F + fc]);
		  //printf("Recovered Uh (t=%d, g=%d, fc=%d) = %f\n",t,g,fc, TestUhCommon[t*M_G*N_F + g*N_F + fc]);
		}
	    }
	}
  
      for (int om = 0; om < N_E; om++)
	{
	  //Start with summation over the element's degrees of freedom.
	  for (int a = 0; a < D; a++)
	    {
	      for (int k = 0; k < N_s; k++)
		{
		  int fc = 0;
		  scalar sum = 0;
		  for (int j = 0; j < N_s; j++)
		    {
		      sum += serial_AuxHatVol[om*D*N_s*N_s + a*N_s*N_s + k*N_s + j] * h_U[(om*N_F+fc)*N_s+j];
		    }
		  TestAuxHat[om*D*N_s*N_F + a*N_s*N_F + k*N_F + fc] = sum;
		}
	    }
	  //Now add the contributions from N_N surrounding elements
	  scalar UhCommon_Element[N_N*M_G*N_F];
	  for (int s = 0; s < N_N; s++) //I think that when true boundary is involved, only need N_N s entries
	    {
	      //int tSide = RecoPair[om*(D+N_N) + s]; //interface we are loading information from
	      int tSide = Alt_FaceFromElem[om*N_N + s]; //interface we are loading information from
	      for (int j = 0; j < M_G*N_F; j++)
		{
		  UhCommon_Element[s*M_G*N_F + j] = TestUhCommon[tSide*M_G*N_F + j];
		}
	    }
	  //Now, we know the UhCommon distribution around the entire element. Use it to 
	  //finish populating GradCommon
	  for (int s = 0; s < N_N; s++) //contributing side of the boundary element
	    {
	      for (int k = 0; k < N_s; k++) //coefficient in auxiliary expansion
		{
		  int fc = 0;
		  for (int a = 0; a < D; a++) //gradient component
		    {
		      scalar sum = 0;
		      for (int j = 0; j < M_G; j++) //contributing quadrature point from face s of omA
			{
			  sum += serial_AuxHatSurf[om*D*N_s*N_N*M_G + a*N_s*N_N*M_G + k*N_N*M_G + s*M_G + j] * UhCommon_Element[s*M_G*N_F + j*N_F + fc];
			}
		      TestAuxHat[om*D*N_s*N_F + a*N_s*N_F + k*N_F + fc] += sum;
		    }
		} //end g loop
	    } //end s loop on faces of boundary element
	}
      printf("pr=%d: The auxiliary coefficients from initial data:\n",myid);
      for (int om = 0; om < N_E; om++)
	{
	  //Start with summation over the element's degrees of freedom.
	  for (int a = 0; a < D; a++)
	    {
	      for (int k = 0; k < N_s; k++)
		{
		  int fc = 0;
		  printf("pr=%d: SigmaHat(e=%d,a=%d,k=%d,fc=%d) = %f\n",myid,om,a,k,fc,TestAuxHat[om*D*N_s*N_F + a*N_s*N_F + k*N_F + fc]);
		}
	    }
	}
    }

  // return 0;


  delete[] serial_AuxHatSurf;
  delete[] serial_AuxHatVol;
  delete[] TestUhCommon;
  delete[] TestAuxHat;
  
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  scalar* TestSigCommon = new scalar[M_T*M_G*N_F*D];
  scalar* TestSigCommonHalf = new scalar[M_T*M_G*N_F*D*2];
  for (int j = 0; j < M_T*M_G*N_F*D*2; j++)
    {
      TestSigCommonHalf[j] = 0.0;
    }
  //test the interface gradient solve
  if (SigFaceTest == 1)
    {
      printf("pr=%d, EXECUTING Face Gradient SOLVE TEST\n",myid);
      for (int e = 0; e < Ne_AUG; e++)
	{
	  scalar sum = 0.0;
	  //scalar FaceRelay[M_G*N_F*D];
	  for (int s = 0; s < D+N_N; s++)
	    {
	      int tGlo = RecoPair[e*(D+N_N) + s]; 
	      int ArgSide = BinarySideAddress[e*(D+N_N) + s]; //either zero or 1, tells the interface which side e is on
	      int fetch_op = e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s); //neighborhood in the SigFace operator
	      int fetch_U = e*N_F*N_s; //neighborhood in the U vector
	      int fetch_out = tGlo*(M_G*N_F*D*2) + ArgSide; //neighborhood in the gradCommonHalf vector
	      //printf("e=%d,s=%d,tGlo=%d,ArgSide=%d\n", e, s, tGlo, ArgSide);
	      int index = 0;
	      for (int g = 0; g < M_G; g++)
		{
		  for(int fc = 0; fc < N_F; fc++)
		    {
		      for (int a = 0; a < D; a++)
			{
			  for (int k = 0; k < N_s; k++)
			    {
			      sum += SigFace_from_Uhat_elemwise[fetch_op + g*(D*N_s) + a*N_s + k] * h_U[fetch_U + fc*N_s + k];
			    }
			  //gotta use += here because many elements send zeros to the zero interface
			  TestSigCommonHalf[fetch_out + g*(N_F*D*2) + fc*D*2 + a*2] += sum;
			  //gradCommonHalf[tGlo*(M_G*N_F*D*2) + g*(N_F*D*2) + fc*D*2 + a*2 + ArgSide] += sum;
			  
			  if (a > 0 && sum > 0.1)
			    {
			      //	      printf("pr=%d:  e=%d/%d, s=%d,g=%d,a=%d,argside=%d); sum to gradCommonHalf(%d,%d,%d,a=%d,side=%d) = %f\n",myid, e,N_E,s,g,a,ArgSide,tGlo,g,fc,a,ArgSide,TestSigCommonHalf[tGlo*(M_G*N_F*D*2) + g*(N_F*D*2) + fc*D*2 + a*2 + ArgSide]);
			    }
			  sum = 0;
			} //end derivative component looop
		    } //end field variable loop
		} //end quadrature node loop
	    }
	}
      //That's half of the gradient on each interface
      for (int t = 0; t < M_T; t++)
	{
	  for (int g = 0; g < M_G; g++)
	    {
	      for (int fc = 0; fc < N_F; fc++)
		{
		  for( int a = 0; a < D; a++)
		    {
		      TestSigCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a] = TestSigCommonHalf[t*M_G*N_F*D*2 + g*N_F*D*2 + fc*D*2 + a*2 + 0] + TestSigCommonHalf[t*M_G*N_F*D*2 + g*N_F*D*2 + fc*D*2 + a*2 + 1];
		      if (fabs( TestSigCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a]) > 0.1 && a > 0)
			{
			  printf("pr=%d: t=%d, (omA=%d/%d, omB=%d/%d), g=%d, a=%d: TestSigCommon=%f\n",myid,t,BR2_Map[t*4+0],N_E, BR2_Map[t*4+2],N_E, g,a,TestSigCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a]);
			}
		      //printf("GradCommon(t=%d,g=%d,f=%d,a=%d) = %f\n", t,g,fc,a,GradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a]);
		    }
		} //end field variable loop
	    } //end quadrature node loop
	} //end interface loop

      //Also, try solving using the interface-specific matrices
      for (int t = 0; t < M_T; t++)
	{
	  int omA = BR2_Map[t*4 + 0];
	  int omB = BR2_Map[t*4 + 2];
	  if (myid==0 && omA == 6 && omB == 32)
	    {
	      printf("SPECIAL INSPECTION at pr=%d, omA=%d, omB=%d\n",myid,omA,omB);
	      for (int k = 0; k < N_s; k++)
		{
		  int fc = 0;
		  printf("UA[index=%d] = %f | ", k, h_U[omA*N_F*N_s + fc*N_s + k]);
		  printf("XYZ_A[index=%d] = (%f, %f, %f)\n", k, XYZNodes_extended[omA*D*N_s + 0*N_s + k], XYZNodes_extended[omA*D*N_s + 1*N_s + k], XYZNodes_extended[omA*D*N_s + 2*N_s + k]);
		}
	      for (int k = 0; k < N_s; k++)
		{
		  int fc = 0;
		  printf("UB[index=%d] = %f | ", k, h_U[omB*N_F*N_s + fc*N_s + k]);
		  printf("XYZ_B[index=%d] = (%f, %f, %f)\n", k, XYZNodes_extended[omB*D*N_s + 0*N_s + k], XYZNodes_extended[omB*D*N_s + 1*N_s + k], XYZNodes_extended[omB*D*N_s + 2*N_s + k]);
		}
	    }
	  for (int g = 0; g < M_G; g++)
	    {
	      for (int a = 0; a < D; a++)
		{
		  scalar sum = 0;
		  int fc = 0;
		  for (int k = 0; k < N_s; k++)
		    {
		    
		      sum += serial_SigFace_from_DOF[t*(D*M_G*2*N_s) + a*M_G*2*N_s + g*2*N_s + k] * h_U[(omA*N_F+fc)*N_s+k];
		      sum += serial_SigFace_from_DOF[t*(D*M_G*2*N_s) + a*M_G*2*N_s + g*2*N_s + N_s + k] * h_U[(omB*N_F+fc)*N_s+k];
		    }
		  printf("pr=%d: Sigface(t=%d, a=%d, g=%d) = %f while TestSigCommon=%f | ", myid, t, a, g, sum, TestSigCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a]);
		}
	      printf("\n");
	    }
	} //end interface (t) loop
    } //end the test for SigFace solve

  //Now safe to delete the serial_SigFace operator
  delete[] serial_SigFace_from_DOF;
  delete[] XYZNodes_extended; //no more use for XYZNodes_extended, it is in matrix form now.


#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif   

  //  printf("pr=%d: About to calculate ADeps values\n",myid);
   
  //Next: Considerations for artificial dissipation.
  //For DGsolver, ADeps only needs to be N_E*D. However,
  //I need to calculate ADepsF for each interface, so ADeps
  //must be calculated for the ghost elements.
  scalar* ADeps = new scalar[Ne_AUG*D]; makeZero(ADeps, Ne_AUG*D); //(h/p) per element
  scalar* ADepsF = new scalar[M_T*D]; makeZero(ADepsF, M_T*D); //(h/p) per interface

  //Now, get the (h/p)^2 parameter per element and also send to APepsF
  for (int e = 0; e < Ne_AUG; e++)
    {
      //This routine copied from setDx in mesh/simpleMesh.cc,
      //except that I want the maximum width per element
#ifdef ONED
      {
	//scalar _Dx = fabs(XYZNodes(1,e) - XYZNodes(0,e)); 
	scalar _Dx = fabs(XYZNodes_extended_Matrix(1, e*D+0) - XYZNodes_extended_Matrix(0, e*D+0));
	//	scalar _Dx = fabs(XYZNodes_extended[e*D*N_s + 0*N_s + 1] - XYZNodes_extended[e*D*N_s + 0*N_s + 0]);
	ADeps[e] = fabs(_Dx / (order+0.0));
      }
#endif
#ifdef TWOD
      {
	//Initialize minimum element width with first edge
	//Parting with Marc's approach here. Given any two possible nodes in the element,
	//I want DelxMax and DelyMax
	scalar DelXMax = 0.0;
	scalar DelYMax = 0.0;
	//For 2D cases, N_N=number of sides = number of faces (not the case in 3D)
	for (int j1 = 0; j1 < N_N; j1++) //hopefully, first N_N nodes are the vertices regardless of quad vs simplex
	  {
	    for (int j2 = 0; j2 < N_N; j2++)
	      {
		DelXMax = fmax(DelXMax, fabs(XYZNodes_extended_Matrix(j1, e*D+0) - XYZNodes_extended_Matrix(j2, e*D+0) ) );
		DelYMax = fmax(DelYMax, fabs(XYZNodes_extended_Matrix(j1, e*D+1) - XYZNodes_extended_Matrix(j2, e*D+1) ) );
		//DelXMax = fmax(DelXMax, fabs(XYZNodes_extended[e*D*N_s + 0*N_s + j1] - XYZNodes_extended[e*D*N_s + 0*N_s + j2]));
		//DelYMax = fmax(DelYMax, fabs(XYZNodes_extended[e*D*N_s + 1*N_s + j1] - XYZNodes_extended[e*D*N_s + 1*N_s + j2]));
		//DelXMax = fmax(DelXMax, fabs(XYZNodes(j1,e*D+0) - XYZNodes(j2,e*D+0)));
		//DelYMax = fmax(DelYMax, fabs(XYZNodes(j1,e*D+1) - XYZNodes(j2,e*D+1)));
	      }
	  }
	ADeps[e*D+0] = DelXMax / (order + 0.0);
	ADeps[e*D+1] = DelYMax / (order + 0.0);
      }
#endif
#ifdef THREED
      {
      // distance from point to surface?
      //Sorry Marc, I don't know either. To start, it's going to be
      //the maximum edge length.
      //Future reference: there is a point to surface distance routine
      //in recoverytools, could bring it in later
      //if I want to be more correct on deformed elements
      if (N_N != 6)
	{
	  printf("CATASTROPHE!!! getDx function not ready for non-hex in 3D\n");
	}
      else
	{
	  scalar DelXMax = 0.0;
	  scalar DelYMax = 0.0;
	  scalar DelZMax = 0.0;
	  //Each hexahedron is defined by 8 vertex nodes; if it is a higher order element, the first eight still form the corners.
	  for (int j1 = 0; j1 < 8; j1++) 
	    {
	      for (int j2 = 0; j2 < 8; j2++)
		{
		  DelXMax = fmax(DelXMax, fabs(XYZNodes_extended_Matrix(j1, e*D+0) - XYZNodes_extended_Matrix(j2, e*D+0) ) );
		  DelYMax = fmax(DelYMax, fabs(XYZNodes_extended_Matrix(j1, e*D+1) - XYZNodes_extended_Matrix(j2, e*D+1) ) );
		  DelZMax = fmax(DelZMax, fabs(XYZNodes_extended_Matrix(j1, e*D+2) - XYZNodes_extended_Matrix(j2, e*D+2) ) );
		  //DelXMax = fmax(DelXMax, fabs(XYZNodes(j1,e*D+0) - XYZNodes(j2,e*D+0)));
		  //DelYMax = fmax(DelYMax, fabs(XYZNodes(j1,e*D+1) - XYZNodes(j2,e*D+1)));
		  //DelZMax = fmax(DelZMax, fabs(XYZNodes(j1,e*D+2) - XYZNodes(j2,e*D+2)));
		}
	    }
	  ADeps[e*D+0] = DelXMax / (order + 0.0);
	  ADeps[e*D+1] = DelYMax / (order + 0.0);
	  ADeps[e*D+2] = DelZMax / (order + 0.0);
	}
      }
#endif
  if (verbose == 1)
    {
      printf("Element %d: ADeps = %f\n", e, ADeps[e]);
    }
} //End element loop to get AD epsilon for each element

  //For interfaces, take Max ADeps from the two nieghboring elements
  for (int t = 0; t < M_T; t++)
    {
      int omA = BR2_Map[t*4 + 0*2 + 0];
      int omB = BR2_Map[t*4 + 1*2 + 0];
      for (int a = 0; a < D; a++)
	{
	  ADepsF[t*D + a] = fmax(ADeps[omA*D + a], ADeps[omB*D + a]);
	   if (verbose == 1)
	     {
	       printf("pr=%d, Interface %d: (omA,omB) = (%d,%d), ADepsF[a=%d] = %f\n",myid, t, omA,omB,a,ADepsF[t*D + a]);
	     }
	}
    }

 

  

   //return 0;
  //
 

  //////////////////////////////////////////////////////////////////////////   
  //
  // Printer setup
  //
  //////////////////////////////////////////////////////////////////////////
  PRINTER printer(N_s,N_E,elem_type,m,timers,mem_counter);
  printer.set_names();

  //////////////////////////////////////////////////////////////////////////   
  //
  // Sensor setup
  //
  //////////////////////////////////////////////////////////////////////////
  SENSOR sensor(N_s, N_E, Ne_AUG, N_N, timers, mem_counter, inputs.getThresholds());

  //////////////////////////////////////////////////////////////////////////   
  //
  // Lagrange particles setup
  //
  //////////////////////////////////////////////////////////////////////////
  LAGRANGE_PARTICLES particles(timers, mem_counter, m, XYZNodes, N_T, N_N, N_E, N_s, myid, numprocs, inputs.getLagrangeParticles());

  //////////////////////////////////////////////////////////////////////////   
  //
  // Setup the dg solver
  //
  //////////////////////////////////////////////////////////////////////////   
  scalar* h_weight  = new scalar[N_G]; makeZero(h_weight,N_G); for(int g=0; g<N_G; g++){h_weight[g] = (scalar)weight(g,0);} mem_counter.addToCPUCounter(N_G*sizeof(scalar));
  printf("processor %d: About to declare dgsolver\n",myid);
  //return 0;
  //BR2 Edit 01/05/2016: Adding extra arguments to the DG_solver call
  DG_SOLVER dgsolver = DG_SOLVER(N_E, N_s, N_G, N_N, M_T, M_s, M_G, M_B, Ne_AUG, communicator, myid, 
  				 h_map, h_invmap, h_phi, h_dphi, h_phi_w, h_dphi_w, h_psi, h_psi_w, //h_xyz, h_xyzf,
  				 h_J, h_invJac, h_JF, h_weight, h_normals, m, timers, mem_counter,
				 /*h_Jac_trace, h_invJac_trace, h_phigF, h_dphigF,*/ BR2_Map, /*h_Minv, h_weight_trace,*/
				 PSIxR_global, PSIxR_elemwise, SigFace_from_Uhat_elemwise, RecoPair, BinarySideAddress, serial_MSigVol, serial_MSigSurf, Alt_FaceFromElem, serial_SigFace_from_Uhat_ROG, PSIxR_biased_global, PSIxR_biased_elemwise, detJxW_full, detJ_full, serial_Uhat2GradBC, serial_UhCommon2GradBC, ADeps, ADepsF, Beta_S, epsGen, Mew_S/*, sensor*/, Cthresh);
  printf("processor %d: Declared dgsolver\n",myid); fflush(stdout);
  
  //return 0;
  //////////////////////////////////////////////////////////////////////////   
  //
  // Setup the limiter
  //
  //////////////////////////////////////////////////////////////////////////
#ifdef ONED
  Limiting Limiter = Limiting(limiterMethod, N_s, N_E, N_N, m, Lag2Mono, Mono2Lag, timers, mem_counter);
#elif TWOD

  Limiting Limiter = Limiting(limiterMethod, N_s, N_E, order, cartesian, N_N, N_G, N_ghosts, m, refArea, Lag2MonoX, MonoX2MonoY, MonoY2Lag, timers, mem_counter);

  // //
  // // Unstructured mesh (relic of the past)
  // //
  // if(!cartesian){
  //   int L = getTaylorDerIdx2DLength(order);
  //   int* h_TaylorDxIdx = new int[L];
  //   int* h_TaylorDyIdx = new int[L];
  //   getTaylorDerIdx2D(order, h_TaylorDxIdx, h_TaylorDyIdx);
  //   Limiting Limiter = Limiting(limiterMethod, N_s, N_E, N_F, N_G, N_N, L, order, L2Msize1, L2Msize2, m.getNeighbors(), Lag2Mono, Mono2Lag, XYZCen, h_powersXYZG, h_weight, refArea, h_TaylorDxIdx, h_TaylorDyIdx);
  //   delete[] h_TaylorDxIdx; delete[] h_TaylorDyIdx;
  //   delete[] h_powersXYZG;
  // }
 
#endif
#ifdef THREED

  Limiting Limiter = Limiting(limiterMethod, N_s, N_E, order, cartesian, N_N, N_G, N_ghosts, m, refArea, Lag2MonoX, MonoX2MonoY, MonoY2Lag, timers, mem_counter);
#endif

  //////////////////////////////////////////////////////////////////////////   
  //
  // Solve the problem on the CPU/GPU. Time integration
  //
  //////////////////////////////////////////////////////////////////////////   
  double DtOut = inputs.getOutputTimeStep();
  double Tf    = inputs.getFinalTime();
  m.setDx(N_N,N_E,XYZCen,XYZNodes);
  //scalar CFL   = inputs.getCFL()*m.getDx()/(2.0*order+1);
  scalar CFL   = inputs.getCFL()*m.getDx();
  scalar VNN = inputs.getVNN()*m.getDx()*m.getDx(); //PEJ 06/01/2017
  scalar VNNAD = inputs.getVNN()*m.getDx()*m.getDx(); //PEJ 10/17/2017
  printf("pr=%d About to Start: input VNN=%f,  VNN w/ geometry=%f, hmin=%f\n", myid, inputs.getVNN(), VNN, m.getDx());
  if (GradMethod_Vis == 0 && HAG == 0)
    {
      switch(order) //The CFL and VNN set
	{
	  //PEJ: no curve fit has yet been calculated for CGR-I (chi=2) timestep restrictions.
	  //So, this will be RK1 values multiplied by appropriate RK4 correction factor.
	  //2D results use shear=1/2 Fourier analysis.
	  //Also making room here for BR2 timestep restriction instead
      
#ifdef ONED
	case 0: VNN = VNN / 4.0;   CFL = CFL / 6.0; break;
	case 1: VNN = VNN / 24.0;  CFL = CFL / 6.0; break;
	case 2: VNN = VNN / 76.5;  CFL = CFL / 12.0; break;
	case 3: VNN = VNN / 179.5; CFL = CFL / 20.0; break;
	case 4: VNN = VNN / 346.5; CFL = CFL / 28.0; break;
	case 5: VNN = VNN / 585.0; CFL = CFL / 38.0; break;
#endif
#ifdef TWOD
	case 0: VNN = VNN / 48.0;   CFL = CFL / 6.0; break;
	case 1: VNN = VNN / 48.0;   CFL = CFL / 6.0; break; 
	case 2: VNN = VNN / 154.0;  CFL = CFL / 12.0; break; 
	case 3: VNN = VNN / 364.0;  CFL = CFL / 20.0; break; 
	case 4: VNN = VNN / 727.0;  CFL = CFL / 28.0; break;
	case 5: VNN = VNN / 1288.0; CFL = CFL / 38.0; break;
#endif
#ifdef THREED
	case 0: VNN = VNN / pow(2, 2) / 4.0;   CFL = CFL / 6.0; break;
	case 1: VNN = VNN / pow(2, 2) / 24.0;  CFL = CFL / 6.0; break; 
	case 2: VNN = VNN / pow(2, 2) / 77.0;  CFL = CFL / 12.0; break;
	case 3: VNN = VNN / pow(2, 2) / 180.0; CFL = CFL / 20.0; break;
	case 4: VNN = VNN / pow(2, 2) / 347.0; CFL = CFL / 28.0; break;
	case 5: VNN = VNN / pow(2, 2) / 585.0; CFL = CFL / 38.0; break; 
#endif
	default: 
	  {
	    VNN = VNN / 5000;
	    printf("\n\nWARNING! p=%d out of range for CGR. Spectral radius unknown\n\n",order);
	    break;
	  }
	}
    }
  else //Either BR2 or HAG, I'll just use BR2 restriction for both of them
    {
      switch(order) //The CFL and VNN set
	{
#ifdef ONED
	case 0: VNN = VNN / 36.0;   CFL = CFL / 6.0; break;
	case 1: VNN = VNN / 36.0;  CFL = CFL / 6.0; break;
	case 2: VNN = VNN / 147.0;  CFL = CFL / 12.0; break;
	case 3: VNN = VNN / 420.0; CFL = CFL / 20.0; break;
	case 4: VNN = VNN / 976.0; CFL = CFL / 28.0; break;
	case 5: VNN = VNN / 1965.0; CFL = CFL / 38.0; break;
#endif
#ifdef TWOD
	case 0: VNN = VNN / 72.0;   CFL = CFL / 6.0; break;
	case 1: VNN = VNN / 72.0;   CFL = CFL / 6.0; break; 
	case 2: VNN = VNN / 294.0;  CFL = CFL / 12.0; break; 
	case 3: VNN = VNN / 841.0;  CFL = CFL / 20.0; break; 
	case 4: VNN = VNN / 1957.0;  CFL = CFL / 28.0; break;
	case 5: VNN = VNN / 3945.0; CFL = CFL / 38.0; break;
#endif
#ifdef THREED
	case 0: VNN = VNN / pow(2, 2) / 36.0;   CFL = CFL / 6.0; break;
	case 1: VNN = VNN / pow(2, 2) / 36.0;  CFL = CFL / 6.0; break; 
	case 2: VNN = VNN / pow(2, 2) / 147.0;  CFL = CFL / 12.0; break;
	case 3: VNN = VNN / pow(2, 2) / 420.0; CFL = CFL / 20.0; break;
	case 4: VNN = VNN / pow(2, 2) / 976.0; CFL = CFL / 28.0; break;
	case 5: VNN = VNN / pow(2, 2) / 1965.0; CFL = CFL / 38.0; break; 
#endif
	default: 
	  {
	    VNN = VNN / 5000;
	    printf("\n\nWARNING! p=%d out of range for CGR. Spectral radius unknown\n\n",order);
	    break;
	  }
	}
    }
  //Account for ICB-N vs. SIMPLE: using 1D analysis from SciTECH draft
#ifdef ICBN
  if (M_B == 0) //only amplify CFL in periodic BC case.
    {
      printf("pr=%d, amplifying CFL for lack of physical boundaries (this is a good thing)\n",myid);
      switch(order)
	{
	case 0: CFL = CFL * 6.0 / 1.3; break;
	case 1: CFL = CFL * 6.0 / 1.3; break;
	case 2: CFL = CFL * 12.0 / 2.7; break;
	case 3: CFL = CFL * 20.0 / 4.2; break;
	case 4: CFL = CFL * 28.0 / 5.8; break;
	case 5: CFL = CFL * 38.0 / 7.4; break;
	default:{ printf("p chosen too high, exiting\n"); exit(1); break;}
	}
    }
#endif
  /*
  switch(order) //THE VNNAD set (CGR-R, Chi=2 values)
    {
    case 0: VNNAD = VNNAD / 4.0; break;
    case 1: VNNAD = VNNAD / 24.0; break;
    case 2: VNNAD = VNNAD / 76.5; break;
    case 3: VNNAD = VNNAD / 179.5; break;
    case 4: VNNAD = VNNAD / 346.5; break;
    case 5: VNNAD = VNNAD / 585.0; break;
    }
  
  //For AD VNN: Account for spatial dimensions. the VNN set
  //accounteed for this within its switch
  VNNAD = VNNAD / pow(2,(D-1));
  */
  VNNAD = VNN; //Presently using same scheme for vis and AD physics, so use same stability restriction

  //Account for RK1 Fourier footprint:
  CFL   = 2.0 * CFL;
  VNN   = 2.0 * VNN; 
  VNNAD = 2.0 * VNNAD;


  if (verbose > 0) {printf("pr=%d, Applied p correction, geometric VNN is now %f\n", myid, VNN);}

  switch (RKchoice)
    {
    case 1:
      {
	RK rk1 = RK(1,1,DtOut,Tf, inputs.getOutputTimeArray());
	printf("pr=%d: ==== Now RK 1 steps with CFL/h = %f, VNN/(h^2) = %f, target time = %e =====\n",myid,CFL/m.getDx(), VNN/(m.getDx()*m.getDx()), Tf);
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	rk1.RK_integration(CFL, VNN, VNNAD, restart_step,
			   N_E, N_s, N_G, M_T, M_s, N_ghosts, N_N, 
			   h_Minv, 
			   h_U, m.getNeighbors(),
			   Limiter, order0, dgsolver, communicator, printer, sensor, timers, mem_counter,particles);
	//printf("pr=%d, RK1 integraion concluded, about to break loop\n",myid);
	break;
      }
      case 2:
      {
        CFL = CFL * 2.0 / 2.0; //RK4 vs RK1 correction 
      VNN = VNN * 2.0 / 2.0; //RK4 vs RK1 correction                                                                                                                         
      VNNAD = VNNAD * 2.0 / 2.0; //RK4 vs RK1 correction                                                                                                                     
        RK rk2 = RK(2,2,DtOut,Tf, inputs.getOutputTimeArray());
        // RK integration                                                                                                                                                    
        printf("pr=%d: ==== Now RK 2 steps with CFL/h = %f, VNN/(h^2) = %f, target time = %e =====\n",myid, CFL/m.getDx(), VNN/(m.getDx()*m.getDx()), Tf);
#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        //PEJ 06/01/2017: addint CNN to the RK4 call                                                                                                                         
        rk2.RK_integration(CFL, VNN, VNNAD, restart_step,
                           N_E, N_s, N_G, M_T, M_s, N_ghosts, N_N,
                           h_Minv,
                           h_U, m.getNeighbors(),
                           Limiter, order0, dgsolver, communicator, printer, sensor, timers, mem_counter,particles);
        //printf("pr=%d, RK4 integraion concluded, about to break loop\n",myid);                                                                                  
        break;
      }
      case 3:
      {
       	CFL = CFL * 2.0 / 2.0; //RK4 vs RK1 correction                                                                                                                       
      VNN = VNN * 2.0 / 2.0; //RK4 vs RK1 correction                                                                                                                         
      VNNAD = VNNAD * 2.0 / 2.0; //RK4 vs RK1 correction                                                                                                                     
        RK rk3 = RK(3,3,DtOut,Tf, inputs.getOutputTimeArray());
        // RK integration                                                                                                                                                    
        printf("pr=%d: ==== Now RK 2 steps with CFL/h = %f, VNN/(h^2) = %f, target time = %e =====\n",myid, CFL/m.getDx(), VNN/(m.getDx()*m.getDx()), Tf);
#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        //PEJ 06/01/2017: addint CNN to the RK4 call                                                                                                                         
        rk3.RK_integration(CFL, VNN, VNNAD, restart_step,
                           N_E, N_s, N_G, M_T, M_s, N_ghosts, N_N,
                           h_Minv,
                           h_U, m.getNeighbors(),
                           Limiter, order0, dgsolver, communicator, printer, sensor, timers, mem_counter,particles);
        //printf("pr=%d, RK4 integraion concluded, about to break loop\n",myid);                                                                                             
        break;
      }
    case 4:
      {
	CFL = CFL * 2.8 / 2.0; //RK4 vs RK1 correction
      VNN = VNN * 2.8 / 2.0; //RK4 vs RK1 correction
      VNNAD = VNNAD * 2.8 / 2.0; //RK4 vs RK1 correction
	RK rk4 = RK(4,4,DtOut,Tf, inputs.getOutputTimeArray());
	// RK integration
        printf("pr=%d: ==== Now RK 4 steps with CFL/h = %f, VNN/(h^2) = %f, target time = %e =====\n",myid, CFL/m.getDx(), VNN/(m.getDx()*m.getDx()), Tf);
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	//PEJ 06/01/2017: addint CNN to the RK4 call
	rk4.RK_integration(CFL, VNN, VNNAD, restart_step,
			   N_E, N_s, N_G, M_T, M_s, N_ghosts, N_N,
			   h_Minv, 
			   h_U, m.getNeighbors(),
			   Limiter, order0, dgsolver, communicator, printer, sensor, timers, mem_counter,particles);
	//printf("pr=%d, RK4 integraion concluded, about to break loop\n",myid);
	break;
      }
    case 8:
      {
	CFL = CFL * 5.0 / 2.0; //RK8 vs RK1 correction
	VNN = VNN * 5.0 / 2.0; //RK8 vs RK1 correction
	VNNAD = VNNAD * 5.0 / 2.0; //RK8 vs RK1 correction
	RK rk8 = RK(8,13,DtOut,Tf, inputs.getOutputTimeArray());
	//return 0;
	// RK integration
	printf("pr=%d:==== Now RK[8/13] steps with CFL/h = %f, VNN/(h^2) = %f, target time = %e =====\n",myid,CFL/m.getDx(), VNN/(m.getDx()*m.getDx()), Tf);
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	//PEJ 06/01/2017: addint CNN to the RK4 call
	rk8.RK_integration(CFL, VNN, VNNAD, restart_step,
			   N_E, N_s, N_G, M_T, M_s, N_ghosts, N_N,
			   h_Minv, 
			   h_U, m.getNeighbors(),
			   Limiter, order0, dgsolver, communicator, printer, sensor, timers, mem_counter,particles);
	//printf("pr=%d, RK8 integraion concluded, about to break loop\n",myid);
	break;
      }
    default:
      {
	printf("pr=%d: RKchoice chosen outside the code's capability. Try again by setting RKchoice to 1,4 or 8 at top of main.cc\n",myid);
	break;
      }
    }

 
  printf("pr=%d: Finished with time integration, now maybe calculate error and close program\n",myid);

  //////////////////////////////////////////////////////////////////////////   
  //
  // Error calcuations
  //
  //////////////////////////////////////////////////////////////////////////

#ifdef ERROR

//build some geometry things
scalar* h_SuperPhi = new scalar[N_E*GQsuper*N_s]; 
scalar* h_SuperJ = new scalar[N_E*GQsuper];
SuperPhi.copyMatrixToPointer(h_SuperPhi);
for (int e = 0; e < N_E; e++){
  //sort of taking a guess here
  for (int g = 0; g < GQsuper; g++)
    {
      h_SuperJ[e*GQsuper + g] = SuperdetJ(e,g); }}

  // Initial condition
  fullMatrix<scalar> Uinit(N_s, N_E*N_F);
#ifdef SCALARAD
  if(sinphil) init_dg_sinphil_singlefluid(N_s, N_E, XYZNodes, Uinit); //PEJ 05/23/2017
  else if(sinphilProject) init_dg_sinphilProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, Uinit); //PEJ 07/05/2017
#endif
#ifdef SINGLEFLUID
  if(tranvtx) init_dg_tranvtx_singlefluid(N_s, N_E, XYZNodes, XYZCen, Uinit, inputs.getInitialConditionInputs());
  else if(HiOWvtx) init_dg_HiOWvtx_singlefluid(N_s, N_E, XYZNodes, XYZCen, Uinit, inputs.getInitialConditionInputs()); //PEJ Edit 01/24/2017
  else if(sodphil) init_dg_sodphil_singlefluid(N_s, N_E, XYZNodes, Uinit); //PEJ Edit 01/24/2017
  else if(normvtx) init_dg_normvtx_singlefluid(N_s, N_E, XYZNodes, Uinit); //PEJ Edit 02/09/2017
  else if (sinphil) init_dg_sinphil_singlefluid(N_s, N_E, XYZNodes, Uinit);
  else if(sinphilProject) init_dg_sinphilProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, Uinit); //PEJ 07/05/2017
  else if(normvtxProject) init_dg_normvtxProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, Uinit); //PEJ 09/08/2017
  else if(HiOWvtxProject) init_dg_HiOWvtxProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, Uinit); //PEJ 10/02/2017
  else if(worsvtxProject) init_dg_worsvtxProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, Uinit); //PEJ 10/18/2017
#elif RADSINGLEFLUID
  if     (sodphil) init_dg_sodphil_singlefluid(N_s, N_E, XYZNodes, Uinit);
  else if(explode) init_dg_explode_singlefluid(N_s, N_E, XYZNodes, Uinit);
  else if(worsvtxProject) init_dg_worsvtxProject_singlefluid(N_s, N_E, GQsuper, SuperXYZ_GQ, SuperWeight, SuperdetJ, SuperPhi, h_Minv, Uinit); //PEJ 10/18/2017
#elif MULTIFLUID
  if     (simplew) init_dg_simplew_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(tranvtx) init_dg_tranvtx_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit, inputs.getInitialConditionInputs());
  else if(sodtube) init_dg_sodtube_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(contact) init_dg_contact_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(rhotact) init_dg_rhotact_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(matfrnt) init_dg_matfrnt_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(sinegam) init_dg_sinegam_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(expogam) init_dg_expogam_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(shckint) init_dg_shckint_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit, inputs.getFinalTime());
  else if(shuoshe) init_dg_shuoshe_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(rarecon) init_dg_rarecon_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(sodcirc) init_dg_sodcirc_multifluid(N_s, N_E, XYZNodes, Uinit);
  else if(rminstb) init_dg_rminstb_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(rmmulti) init_dg_rmmulti_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(rtaylor) init_dg_rtaylor_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(khdrake) init_dg_khdrake_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit, inputs.getInitialConditionInputs());
  else if(khuramp) init_dg_khuramp_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit, inputs.getInitialConditionInputs());
  else if(khblast) init_dg_khblast_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(khpertu) init_dg_khpertu_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);
  else if(blastrm) init_dg_blastrm_multifluid(N_s, N_E, XYZNodes, XYZCen, Uinit);

#elif PASSIVE
  if (sinephi) init_dg_sinephi_passive(N_s, N_E, XYZNodes, Uinit);
  if (sodmono) init_dg_sodmono_passive(N_s, N_E, XYZNodes, Uinit);
#endif




  scalar kwave = 2.0*M_PI/DomainLength;

  //decay depends on viscosity:
#ifdef CONSTANTVIS
  scalar mew = constants::GLOBAL_KLIN;
#endif
#ifdef NOVIS
  scalar mew = 0.0;
#endif

  scalar Tau_decay;
#ifdef ONED
  Tau_decay = kwave*kwave*Tf*mew;
#endif
#ifdef TWOD
  Tau_decay = 2.0*kwave*kwave*Tf*mew;
#endif
#ifdef THREED
  Tau_decay = 3.0*kwave*kwave*Tf*mew;
#endif
  printf("Station A: Applying decay factor to initial condition based on viscosity: DomainLength=%f. kwave=%f, mew=%f, Tau_decay=%f\n",DomainLength, kwave, mew, Tau_decay);



  scalar* h_Uinit = new scalar[N_s*N_E*N_F];  makeZero(h_Uinit,N_s*N_E*N_F);
  if(sinphilProject || sinphil)
    {
      for(int e = 0; e < N_E; e++){
      for(int fc = 0; fc < N_F; fc++){
      for(int i = 0; i < N_s; i++){
      //	h_Uinit[(e*N_F+fc)*N_s+i] = Uinit(i,e*N_F+fc);
      //This formula is particular to initial condition in sinphil:
      h_Uinit[(e*N_F+fc)*N_s+i] = 2.0 + exp(-Tau_decay)*(Uinit(i,e*N_F+fc)-2.0);
    }}}
    }
  if (normvtxProject || normvtx || HiOWvtx || HiOWvtxProject)
    {
      for(int e = 0; e < N_E; e++){
	for(int fc = 0; fc < N_F; fc++){
	  for(int i = 0; i < N_s; i++){
	    h_Uinit[(e*N_F+fc)*N_s+i] = Uinit(i,e*N_F+fc);
      }}}
    }
    printf("Station B\n");
    // Change to primitive variables.
  //PEJ 07/05/2017: Commenting out the transformation to primitive variables because
  //it doesn't work for scalar AD case and I don't like it anyway.
  /*
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
      // velocity
      h_Uinit[(e*N_F+1)*N_s+i] = h_Uinit[(e*N_F+1)*N_s+i]/h_Uinit[(e*N_F+0)*N_s+i]; 
      h_U    [(e*N_F+1)*N_s+i] = h_U    [(e*N_F+1)*N_s+i]/h_U    [(e*N_F+0)*N_s+i];
#ifdef SINGLEFLUID
      scalar gamma = constants::GLOBAL_GAMMA;
      h_Uinit[(e*N_F+2)*N_s+i] = (gamma-1)*(h_Uinit[(e*N_F+2)*N_s+i] - 0.5*h_Uinit[(e*N_F+0)*N_s+i]*h_Uinit[(e*N_F+1)*N_s+i]*h_Uinit[(e*N_F+1)*N_s+i]); 
      h_U    [(e*N_F+2)*N_s+i] = (gamma-1)*(h_U    [(e*N_F+2)*N_s+i] - 0.5*h_U    [(e*N_F+0)*N_s+i]*h_U    [(e*N_F+1)*N_s+i]*h_U    [(e*N_F+1)*N_s+i]);
#elif MULTIFLUID
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
*/

    
  // Collocate the solution to the integration points
/*
  scalar* h_Uinitg = new scalar[N_G*N_E*N_F];  makeZero(h_Uinitg,N_G*N_E*N_F);
  scalar* h_Ug     = new scalar[N_G*N_E*N_F];  makeZero(h_Ug    ,N_G*N_E*N_F);
  hostblasGemm('N','N', N_G, N_E*N_F, N_s, 1, h_phi, N_G, h_Uinit, N_s, 0.0, h_Uinitg, N_G);
  hostblasGemm('N','N', N_G, N_E*N_F, N_s, 1, h_phi, N_G, h_U    , N_s, 0.0, h_Ug    , N_G);
*/
  scalar* h_Uinitg = new scalar[GQsuper*N_E*N_F];  makeZero(h_Uinitg,GQsuper*N_E*N_F);
  scalar* h_Ug     = new scalar[GQsuper*N_E*N_F];  makeZero(h_Ug    ,GQsuper*N_E*N_F);
  hostblasGemm('N','N', GQsuper, N_E*N_F, N_s, 1, h_SuperPhi, GQsuper, h_Uinit, N_s, 0.0, h_Uinitg, GQsuper);
  hostblasGemm('N','N', GQsuper, N_E*N_F, N_s, 1, h_SuperPhi, GQsuper, h_U    , N_s, 0.0, h_Ug    , GQsuper);
/*
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
*/
// Take the cell average of the solution
  scalar* h_UinitAvg = new scalar[N_E*N_F];  makeZero(h_UinitAvg,N_E*N_F);
  scalar* h_UAvg     = new scalar[N_E*N_F];  makeZero(h_UAvg    ,N_E*N_F);
  scalar dx = XYZNodes(1,0*D+0)-XYZNodes(0,0*D+0);
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      for(int g = 0; g < GQsuper; g++){
	h_UinitAvg[e*N_F+fc] += h_Uinitg[(e*N_F+fc)*GQsuper+g]*h_SuperJ[e*GQsuper + g]*SuperWeight(g,0);
	h_UAvg    [e*N_F+fc] += h_Ug    [(e*N_F+fc)*GQsuper+g]*h_SuperJ[e*GQsuper + g]*SuperWeight(g,0);
      }
    }
  }
  
  // Calculate the different cell-average norms of the error, Same regardless of quadrature order
  scalar E = 0;
  scalar EMod = 0;
  scalar* h_Err1      = new scalar[N_F]; makeZero(h_Err1   , N_F);
  scalar* h_Err2      = new scalar[N_F]; makeZero(h_Err2   , N_F);
  scalar* h_ErrInf    = new scalar[N_F]; makeZero(h_ErrInf   , N_F);
  for(int e = 0; e < N_E; e++){
    for(int fc = 0; fc < N_F; fc++){
      E = h_UinitAvg[e*N_F+fc]/dx-h_UAvg[e*N_F+fc]/dx;
      h_Err1[fc] += fabs(E);
      h_Err2[fc] += E*E;
      if (h_ErrInf[fc] < fabs(E))  h_ErrInf[fc] = fabs(E);
    }
  }
  printf("Station C\n");
  /*
  //PEJ 07/05/2017: We're going to calculate some more error norms
  scalar* E2G = new scalar[N_F]; makeZero(E2G, N_F);
  scalar* E2J = new scalar[N_F]; makeZero(E2J, N_F);
  scalar* UvsJ_exact = new scalar[N_F]; makeZero(UvsJ_exact, N_F);
  scalar* UvsJ_h = new scalar[N_F]; makeZero(UvsJ_h, N_F);
  //Calculate finite-dimensional interpolation of the exact solution.
  //Involves multiplying IC by diffusive decay term
  
  scalar* h_Uexact = new scalar[N_E*N_s*N_F]; makeZero(h_Uexact, N_E*N_s*N_F);
  scalar* ErrorFunctional = new scalar[N_E*N_s*N_F]; makeZero(ErrorFunctional, N_E*N_s*N_F);
  //populate each Lagrange (solution) node with the exact solution at end of simulation
  for (int e = 0; e < N_E; e++)
    {
      for (int fc = 0; fc < N_F; fc++)
	{
	  for (int k = 0; k < N_s; k++)
	    {
	      h_Uexact[e*N_F*N_s + fc*N_s + k] = 2.0 + exp(-Tau_decay)*sin(2.0*M_PI/10.0*XYZNodes(k,e*D + 0)); //sin(x)
	      ErrorFunctional[e*N_F*N_s + fc*N_s + k] = sin(1.1*2.0*M_PI/10.0*XYZNodes(k,e*D+0));
	    }
	}
    } 
  //send the exact solution and error functional to integration points
  scalar* h_Uexactg = new scalar[N_G*N_E*N_F];  makeZero(h_Uexactg,N_G*N_E*N_F);
  hostblasGemm('N','N', N_G, N_E*N_F, N_s, 1, h_phi, N_G, h_Uexact, N_s, 0.0, h_Uexactg, N_G);
  scalar* ErrorFunctional_g = new scalar[N_E*N_G*N_F]; makeZero(ErrorFunctional_g, N_E*N_G);
  hostblasGemm('N','N', N_G, N_E*N_F, N_s, 1, h_phi, N_G, ErrorFunctional, N_s, 0.0, ErrorFunctional_g, N_G); 

  //Now compare h_Uexactg and h_Ug with some error norms
  for (int e = 0; e < N_E; e++)
    {
      for (int fc = 0; fc < N_F; fc++)
	{
	  for (int g = 0; g < N_G; g++)
	    {
	      E2G[fc] += h_J[e]*weight(g,0)*pow(h_Uexactg[e*N_F*N_G + fc*N_G + g] - h_Ug[e*N_F*N_G + fc*N_G + g],2);
	      UvsJ_exact[fc] += h_J[e]*weight(g,0)*ErrorFunctional_g[e*N_F*N_G + fc*N_G + g]*h_Uexactg[e*N_F*N_G + fc*N_G + g];
	      UvsJ_h[fc] += h_J[e]*weight(g,0)*ErrorFunctional_g[e*N_F*N_G + fc*N_G + g]*h_Ug[e*N_F*N_G + fc*N_G + g];
	    }
	}
    }
  for (int fc = 0; fc < N_F; fc++)
    {
      E2G[fc] = sqrt(E2G[fc]);
      E2J[fc] = fabs(UvsJ_exact[fc] - UvsJ_h[fc]);
    }
  */
  //PEJ 07/05/2017: We're going to calculate some more error norms
  scalar* E2G = new scalar[N_F]; makeZero(E2G, N_F);
  scalar* E2J = new scalar[N_F]; makeZero(E2J, N_F);
  scalar* UvsJ_exact = new scalar[N_F]; makeZero(UvsJ_exact, N_F);
  scalar* UvsJ_h = new scalar[N_F]; makeZero(UvsJ_h, N_F);
  //Calculate finite-dimensional interpolation of the exact solution.
  //Involves multiplying IC by diffusive decay term
  
  scalar* h_Uexact = new scalar[N_E*N_s*N_F]; makeZero(h_Uexact, N_E*N_s*N_F);
  scalar* ErrorFunctional = new scalar[N_E*N_s*N_F]; makeZero(ErrorFunctional, N_E*N_s*N_F);
  //populate each Lagrange (solution) node with the exact solution at end of simulation
  for (int e = 0; e < N_E; e++)
    {
      for (int fc = 0; fc < N_F; fc++)
	{
	  for (int k = 0; k < N_s; k++)
	    {
	      h_Uexact[e*N_F*N_s + fc*N_s + k] = h_Uinit[(e*N_F+fc)*N_s+k]; 
	      //ErrorFunctional[e*N_F*N_s + fc*N_s + k] = sin(1.1*2.0*M_PI/10.0*XYZNodes(k,e*D+0));
	      ErrorFunctional[e*N_F*N_s + fc*N_s + k] = sin(1.1*kwave*XYZNodes(k,e*D+0));
	    }
	}
    } 
  printf("Station D\n");
  //send the exact solution and error functional to integration points
  scalar* h_Uexactg = new scalar[GQsuper*N_E*N_F];  makeZero(h_Uexactg,GQsuper*N_E*N_F);
  //hostblasGemm('N','N', GQsuper, N_E*N_F, N_s, 1, h_SuperPhi, GQsuper, h_Uexact, N_s, 0.0, h_Uexactg, GQsuper); (DO NOT DELETE~ This approach could be useful in future, as it populates not the exact solution, but the finite-dimension projection of the exact solution)
  for (int e = 0; e < N_E; e++)
    {
      for (int g = 0; g < GQsuper; g++)
	{
	  //Get local super-quadrature node location
	  scalar* XYZlocal = new scalar[D];
	  for (int a = 0; a < D; a++)
	    {
	      XYZlocal[a] = SuperXYZ_GQ(g, e*D+a);
	    }
	  //Get exact solution at the quadrature node
	  scalar* sol_local = new scalar[N_F];
	  //Start by setting sol_local=0, then adjust
	  //it if analystical solution is available
	  Exact_zero(XYZlocal, sol_local);
	  if(sinphil)        
	    {
	      Exact_sinphil(XYZlocal, sol_local, DomainLength);
	      sol_local[0] = 2.0 + (sol_local[0]-2.0)*exp(-Tau_decay);
	    }
	  if(sinphilProject) 
	    {
	      Exact_sinphil(XYZlocal, sol_local, DomainLength);
	      sol_local[0] = 2.0 + (sol_local[0]-2.0)*exp(-Tau_decay);
	    }
	  if(normvtx)        Exact_normvtx(XYZlocal, sol_local);
	  if(normvtxProject) Exact_normvtx(XYZlocal, sol_local);
	  if(HiOWvtx) Exact_HiOWvtx(XYZlocal, sol_local);
	  if(HiOWvtxProject) Exact_HiOWvtx(XYZlocal, sol_local);
	  //Store the exact solution at the quadrature node:
	  for (int fc = 0; fc < N_F; fc++)
	    {
	      h_Uexactg[e*N_F*GQsuper + fc*GQsuper + g] = sol_local[fc]; 
	    }
	  delete[] XYZlocal;
	  delete[] sol_local;
	}
    }
  
  printf("Station E\n");
  scalar* ErrorFunctional_g = new scalar[N_E*N_F*GQsuper]; makeZero(ErrorFunctional_g,N_E*N_F*GQsuper);

  //This ErrorFunctional distribution is the one that matters, earlier one is irrelevant
  for (int e = 0; e < N_E; e++) {
    for (int g = 0; g < GQsuper; g++) {
      for (int fc = 0; fc < N_F; fc++) {
#ifdef ONED
	ErrorFunctional_g[e*N_F*GQsuper + fc*GQsuper + g] = sin(1.1*kwave*SuperXYZ_GQ(g,e*D+0)); 
#endif
#ifdef TWOD
	//ErrorFunctional_g[e*N_F*GQsuper + fc*GQsuper + g] = sin(0.25*2.0*M_PI/10.0*SuperXYZ_GQ(g,e*D+0)) + sin(0.75*2.0*M_PI/10.0*SuperXYZ_GQ(g,e*D+1)) + sin(0.5*2.0*M_PI/10.0*SuperXYZ_GQ(g,e*D+0)*SuperXYZ_GQ(g,e*D+1)); 
	ErrorFunctional_g[e*N_F*GQsuper + fc*GQsuper + g] = sin(0.25*kwave*SuperXYZ_GQ(g,e*D+0)) + sin(0.75*kwave*SuperXYZ_GQ(g,e*D+1)) + sin(0.5*kwave*SuperXYZ_GQ(g,e*D+0)*SuperXYZ_GQ(g,e*D+1)); 
#endif
#ifdef THREED
	ErrorFunctional_g[e*N_F*GQsuper + fc*GQsuper + g] = sin(0.25*kwave*SuperXYZ_GQ(g,e*D+0)) + sin(0.75*kwave*SuperXYZ_GQ(g,e*D+1)) + sin(0.5*kwave*SuperXYZ_GQ(g,e*D+0)*SuperXYZ_GQ(g,e*D+1)); 
	//ErrorFunctional_g[e*N_F*GQsuper + fc*GQsuper + g] = sin(0.25*2.0*M_PI/10.0*SuperXYZ_GQ(g,e*D+0)) + sin(1.75*2.0*M_PI/10.0*SuperXYZ_GQ(g,e*D+1)); 
#endif
      }}}

  //Now compare h_Uexactg and h_Ug with some error norms
  for (int e = 0; e < N_E; e++) {
      for (int fc = 0; fc < N_F; fc++) {
	  for (int g = 0; g < GQsuper; g++) {
	    E2G[fc] += h_SuperJ[e*GQsuper + g]*SuperWeight(g,0)*pow(h_Uexactg[e*N_F*GQsuper + fc*GQsuper + g] - h_Ug[e*N_F*GQsuper + fc*GQsuper + g],2);
	    UvsJ_exact[fc] += h_SuperJ[e*GQsuper + g]*SuperWeight(g,0)*ErrorFunctional_g[e*N_F*GQsuper + fc*GQsuper + g]*h_Uexactg[e*N_F*GQsuper + fc*GQsuper + g];
	    UvsJ_h[fc] += h_SuperJ[e*GQsuper + g]*SuperWeight(g,0)*ErrorFunctional_g[e*N_F*GQsuper + fc*GQsuper + g]*h_Ug[e*N_F*GQsuper + fc*GQsuper + g];
	    }
	}
    }
  for (int fc = 0; fc < N_F; fc++)
    {
      E2G[fc] = sqrt(E2G[fc]);
      E2J[fc] = fabs(UvsJ_exact[fc] - UvsJ_h[fc]);
    }
    printf("Station F\n");
  //Also, get cell-average error. h_Err2 is insufficient for my needs
  scalar* E2bar = new scalar[N_F]; makeZero(E2bar, N_F);
  scalar* E2barPrim = new scalar[N_F]; makeZero(E2barPrim, N_F); //cell-average error, primitive variables
  scalar* E2GPrim = new scalar[N_F]; makeZero(E2GPrim, N_F); //global error, primitive variables
  for (int e = 0; e < N_E; e++)
    {
      scalar CellArea = 0.0; //Area of the actual element
      for (int g = 0; g < GQsuper; g++)
	{
	  CellArea += h_SuperJ[e*GQsuper + g] *  SuperWeight(g,0);
	}
      //printf("e=%d, CellArea=%f\n",e,CellArea);
      for (int fc = 0; fc < N_F; fc++)
	{
	  scalar Ubarh = 0;
	  scalar Ubarex = 0;
	  for (int g = 0; g < GQsuper; g++)
	    {
	      Ubarex += h_Uexactg[e*N_F*GQsuper + fc*GQsuper + g] * h_SuperJ[e*GQsuper + g] * SuperWeight(g,0);
	      Ubarh += h_Ug[e*N_F*GQsuper + fc*GQsuper + g] * h_SuperJ[e*GQsuper + g] * SuperWeight(g,0);
	    }
	  Ubarh = Ubarh / CellArea;
	  Ubarex = Ubarex / CellArea;
	  E2bar[fc] += pow(Ubarh - Ubarex,2);
	}
      //Rinse and repeat for primitive variables:
      
#ifdef SINGLEFLUID
      if (D == 2)
	{
	  //Must first caluclate primitive variable distribution at quadrature points
	  scalar rhoG_ex[GQsuper];
	  scalar uG_ex[GQsuper];
	  scalar vG_ex[GQsuper];
	  scalar pG_ex[GQsuper];
	  scalar rhoG_h[GQsuper];
	  scalar uG_h[GQsuper];
	  scalar vG_h[GQsuper];
	  scalar pG_h[GQsuper];
	  for (int g = 0; g < GQsuper; g++)
	    {
	      rhoG_ex[g] = h_Uexactg[e*N_F*GQsuper + 0*GQsuper + g];
	      uG_ex[g]   = h_Uexactg[e*N_F*GQsuper + 1*GQsuper + g] / rhoG_ex[g];
	      vG_ex[g]   = h_Uexactg[e*N_F*GQsuper + 2*GQsuper + g] / rhoG_ex[g];
	      pG_ex[g]   = (1.4-1.0)*(h_Uexactg[e*N_F*GQsuper + 3*GQsuper + g]-0.5*rhoG_ex[g]*(uG_ex[g]*uG_ex[g] + vG_ex[g]*vG_ex[g]));

	      rhoG_h[g] = h_Ug[e*N_F*GQsuper + 0*GQsuper + g];
	      uG_h[g]   = h_Ug[e*N_F*GQsuper + 1*GQsuper + g] / rhoG_h[g];
	      vG_h[g]   = h_Ug[e*N_F*GQsuper + 2*GQsuper + g] / rhoG_h[g];
	      pG_h[g]   = (1.4-1.0)*(h_Ug[e*N_F*GQsuper + 3*GQsuper + g]-0.5*rhoG_h[g]*(uG_h[g]*uG_h[g] + vG_h[g]*vG_h[g]));
	    }
	  //Primitives known at quadrature points. Now, get cell-average error contribution in each primitive field.
	  //Also, while I'm here, might as well get the global L2 error in the primitive variables
	  scalar rhobar_h = 0;
	  scalar rhobar_ex = 0;
	  scalar ubar_h = 0;
	  scalar ubar_ex = 0;
	  scalar vbar_h = 0;
	  scalar vbar_ex = 0;
	  scalar pbar_h = 0;
	  scalar pbar_ex = 0;
	  for (int g = 0; g < GQsuper; g++)
	    {
	      rhobar_ex += rhoG_ex[g] * h_SuperJ[e*GQsuper + g] * SuperWeight(g,0);
	      rhobar_h  += rhoG_h[g]  * h_SuperJ[e*GQsuper + g] * SuperWeight(g,0);

	      ubar_ex   += uG_ex[g]   * h_SuperJ[e*GQsuper + g] * SuperWeight(g,0);
	      ubar_h    += uG_h[g]    * h_SuperJ[e*GQsuper + g] * SuperWeight(g,0);

	      vbar_ex   += vG_ex[g]   * h_SuperJ[e*GQsuper + g] * SuperWeight(g,0);
	      vbar_h    += vG_h[g]    * h_SuperJ[e*GQsuper + g] * SuperWeight(g,0);

	      pbar_ex   += pG_ex[g]   * h_SuperJ[e*GQsuper + g] * SuperWeight(g,0);
	      pbar_h    += pG_h[g]    * h_SuperJ[e*GQsuper + g] * SuperWeight(g,0);

	      E2GPrim[0] += pow(rhoG_ex[g] - rhoG_h[g],2) * SuperWeight(g,0) * h_SuperJ[e*GQsuper + g];
	      E2GPrim[1] += pow(uG_ex[g]   - uG_h[g],2)   * SuperWeight(g,0) * h_SuperJ[e*GQsuper + g];
	      E2GPrim[2] += pow(vG_ex[g]   - vG_h[g],2)   * SuperWeight(g,0) * h_SuperJ[e*GQsuper + g];
	      E2GPrim[3] += pow(pG_ex[g]   - pG_h[g],2)   * SuperWeight(g,0) * h_SuperJ[e*GQsuper + g];
	    }
	  rhobar_h  = rhobar_h / CellArea;
	  rhobar_ex = rhobar_ex / CellArea;

	  ubar_h  = ubar_h / CellArea;
	  ubar_ex = ubar_ex / CellArea;

	  vbar_h  = vbar_h / CellArea;
	  vbar_ex = vbar_ex / CellArea;

	  pbar_h  = pbar_h / CellArea;
	  pbar_ex = pbar_ex / CellArea;

	  //Applying the HiOCFD5 cell-average error norm here: multiply the square of differences by the local element area
	  E2barPrim[0] += pow(rhobar_h - rhobar_ex,2) * CellArea;
	  E2barPrim[1] += pow(ubar_h - ubar_ex,2) * CellArea;
	  E2barPrim[2] += pow(vbar_h - vbar_ex,2) * CellArea;
	  E2barPrim[3] += pow(pbar_h - pbar_ex,2) * CellArea;
	}
    
#endif
    }
  for (int fc = 0; fc < N_F; fc++)
    {
      E2bar[fc] = sqrt(E2bar[fc]/N_E);
      //E2barPrim[fc] = sqrt(E2bar[fc]/N_E);
      E2barPrim[fc] = sqrt(E2barPrim[fc] / pow(DomainLength,D)); //HiOCFD5 definition of cell-average error
      E2GPrim[fc] = sqrt(E2GPrim[fc] / pow(DomainLength,D)); //need to divide global error norm by domain area, then take square root.
    }
  

  // Output some stuff in a file to read by post-proc
  std::string error = "error.dat"; 
  FILE *f = fopen(error.c_str(),"w");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", h_Err1[fc]/(double)N_E);          fprintf(f,"\n");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", sqrt(h_Err2[fc]/(double)N_E));    fprintf(f,"\n");
  fprintf(f,"%12.7f\t", dx); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", h_ErrInf[fc]);                    fprintf(f,"\n");
  fclose(f);

  //Phil's error quantification section
  FILE*file_error = fopen("ErrorPhil.csv","a");
  fprintf(file_error, 
	  "%d, %d,  %d, %f,       %16.13f,",
	  N_E, N_s, D, m.getDx(), pow(N_E*N_s, -1.0/D));
  for (int fc = 0; fc < N_F; fc++)
    {
      fprintf(file_error, "%16.13f,%16.13f,%16.13f,%16.13f,%16.13f,",E2G[fc] , E2bar[fc], E2J[fc], UvsJ_exact[fc], UvsJ_h[fc]);
      //fprintf(file_error, "%16.13f,%16.13f,%16.13f,%16.13f,%16.13f,",E2G[fc] , sqrt(h_Err2[fc]/(double)N_E), E2J[fc], UvsJ_exact[fc], UvsJ_h[fc]);
    }
  fprintf(file_error, "\n");
  fclose(file_error);

#ifdef SINGLEFLUID
  if (D == 2)
    {
      //For the HiOCFD5 vortex case: I've been informed
      //that the quantity of interest is the global L2 errors
      //in u and v.
      //Also need timer information.
      
      FILE*Efile = fopen("ErrorHiOWvort.csv","a");
      fprintf(Efile, 
	  "%d, %d,  %d, %f,       %16.13f,",
	  N_E, N_s, D, m.getDx(), pow(N_E*N_s, -1.0/D));
      for (int fc = 0; fc < N_F; fc++)
	{
	  fprintf(Efile, "%16.13f,%16.13f,",E2barPrim[fc],E2GPrim[fc]);
	}
      //     fprintf(Efile, "%16.8f",timers._times[1]);
      fprintf(Efile, "%16.8f,",timers.print_single_timer(1)); //print total time spent in forward time integration
      fprintf(Efile, "%16.13f,", (timers.print_single_timer(41)+timers.print_single_timer(42))/timers.print_single_timer(1)); //percent of time devoted to ICB routine
      fprintf(Efile, "\n");
      fclose(Efile);
    }
#endif
  
  // Free some stuff
  delete[] h_Uinit; printf("deleted h_Unit\n"); fflush(stdout);
  delete[] h_Uinitg; printf("deleted h_Unitg\n"); fflush(stdout);
  delete[] h_UinitAvg; printf("deleted h_UnitAvg\n"); fflush(stdout);
  delete[] h_Ug; printf("deleted h_Ug\n"); fflush(stdout);
  delete[] h_UAvg; printf("deleted h_UAvg\n"); fflush(stdout);
  delete[] h_Err1; printf("deleted h_Err1\n"); fflush(stdout);
  delete[] h_Err2; printf("deleted h_Err2\n"); fflush(stdout);
  delete[] h_ErrInf; printf("deleted h_ErrInf\n"); fflush(stdout);

  //Phil's additions
  delete[] UvsJ_exact; printf("deleted UvsJ_exact\n"); fflush(stdout);
  delete[] UvsJ_h; printf("deleted UvsJ_h\n"); fflush(stdout);
  delete[] E2J; printf("deleted E2J\n"); fflush(stdout);
  delete[] E2G; printf("deleted E2G\n"); fflush(stdout);
  delete[] h_Uexact; printf("deleted h_Uexact\n"); fflush(stdout);
  delete[] h_Uexactg; printf("deleted h_Uexactg\n"); fflush(stdout);
  delete[] ErrorFunctional; printf("deleted ErrorFunctional\n"); fflush(stdout);
  delete[] ErrorFunctional_g; printf("deleted ErrorFunctional_g\n"); fflush(stdout);
  delete[] h_SuperJ; printf("deleted h_SuperJ\n"); fflush(stdout);
  delete[] h_SuperPhi; printf("deleted h_SuperPhi\n"); fflush(stdout);
  delete[] E2bar; printf("deleted E2bar\n"); fflush(stdout);
  delete[] E2barPrim; printf("deleted E2barPrim\n"); fflush(stdout);
  delete[] E2GPrim; printf("deleted E2GPrim\n"); fflush(stdout);
#endif

  
  //////////////////////////////////////////////////////////////////////////   
  //
  // Free stuff on the device
  //
  //////////////////////////////////////////////////////////////////////////   
#ifdef USE_GPU
  cudaFreeHost(h_U);
  status = cublasShutdown();
#endif
  
  
  //////////////////////////////////////////////////////////////////////////   
  //
  // Free stuff on the host
  //
  //////////////////////////////////////////////////////////////////////////   
delete[] h_Minv; //printf("deleted h_Minv\n"); fflush(stdout);
delete[] h_map; //printf("deleted h_map\n"); fflush(stdout);
  //delete[] h_invmap; printf("deleted h_invmap\n"); fflush(stdout);
delete[] h_phi; //printf("deleted h_phi\n"); fflush(stdout);
delete[] h_phi_w; //printf("deleted h_phi_w\n"); fflush(stdout);
delete[] h_dphi; //printf("deleted h_dphi\n"); fflush(stdout);
delete[] h_dphi_w; //printf("deleted h_dphi_w\n"); fflush(stdout);
delete[] h_psi; //printf("deleted h_psi\n"); fflush(stdout);
delete[] h_psi_w; //printf("deleted h_psi_w\n"); fflush(stdout);
 
  //delete[] h_xyz;
  //delete[] h_xyzf;
delete[] h_weight; //printf("deleted h_weight\n"); fflush(stdout);
delete[] h_J; //printf("deleted h_J\n"); fflush(stdout);
delete[] h_JF; //printf("deleted h_JF\n"); fflush(stdout);
delete[] h_invJac; //printf("deleted h_invJac\n"); fflush(stdout);
delete[] h_normals; //printf("deleted h_normals\n"); fflush(stdout);
 
  //10/24/2017: Looking to fight some segfaults, introducing more deletes
delete[] h_phigF; //printf("deleted h_phigF\n"); fflush(stdout);
delete[] h_dphigF; //printf("deleted h_dphigF\n"); fflush(stdout);
delete[] h_invJac_trace; //printf("deleted h_invJac_trace\n"); fflush(stdout);
delete[] h_Jac_trace; //printf("deleted h_Jac_trace\n"); fflush(stdout);
delete[] h_weight_trace; //printf("deleted h_weight_trace\n"); fflush(stdout);
delete[] FaceFromElem; //printf("deleted FaceFromElem\n"); fflush(stdout);

  //05/27/2017: I guess since Marc deletes stuff here, I should too
delete[] serial_MSigSurf; //printf("deleted serial_MSigSurf\n"); fflush(stdout);
delete[] serial_MSigVol; //printf("deleted seiral_MSigVol\n"); fflush(stdout);
delete[] BR2_Map; //printf("deleted BR2_Map\n"); fflush(stdout);
  //delete[] PSIxR_global;
delete[] PSIxR_elemwise; //printf("deleted PSIxR_elemwise\n"); fflush(stdout);
  


//printf("delete[] PSIxR_biased_elemwise\n"); fflush(stdout);
  delete[] PSIxR_biased_elemwise;
  delete[] Better_InvMap;
//printf("delete[] BinarySideAddress\n"); fflush(stdout);
  delete[] BinarySideAddress;
  delete[] RecoPair;
  delete[] SigFace_from_Uhat_elemwise;
//printf("delete[] Alt_FaceFromElem\n"); fflush(stdout);
  delete[] Alt_FaceFromElem;
  //delete[] serial_SigFace_from_DOF;
  //delete[] UicbA_test;
  //delete[] UicbB_test;
  //delete[] Uicb_test;
  //delete[] UicbHalf_test;
  //delete[] Uhat_fake;
  delete[] detJxW_full;
  delete[] detJ_full;
//printf("delete[] serial_Uhat2GradBC\n"); fflush(stdout);
  delete[] serial_Uhat2GradBC;
//printf("delete[] serial_UhCommon2GradBC\n"); fflush(stdout);
  delete[] serial_UhCommon2GradBC;
//printf("delete[] ADeps\n"); fflush(stdout);
  delete[] ADeps;
  //printf("delete[] ADepF\n"); fflush(stdout);
  delete[] ADepsF;
  delete[] serial_SigFace_from_Uhat_ROG;
  delete[] PSIxR_global;
  delete[] PSIxR_biased_global;
delete[] h_invmap; //printf("deleted h_invmap\n"); fflush(stdout); //causing problems

printf("pr=%d: Last pointer deletion: delete[] h_invmap\n",myid);

#ifdef USE_CPU
  delete[] h_U;
#endif 

  //////////////////////////////////////////////////////////////////////////   
  //
  // End timer and output counters
  //
  //////////////////////////////////////////////////////////////////////////
  timers.stop_timer(0);
  timers.print_timers();
  mem_counter.outputCounters();
  
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
  // And we are done.
  //
  //////////////////////////////////////////////////////////////////////////
  if(myid==0){
    printf("\n\n\n");
    printf("+--------------------------------------------------------------------------------+\n");
    printf("|+------------------------------------------------------------------------------+|\n");
    printf("||                                                                              ||\n");
    printf("||                                                                              ||\n");
    printf("||                                                                              ||\n");
    printf("||                            CALCULATIONS CORRECT.                             ||\n"); 
    printf("||                                                                              ||\n");
    printf("||                                                                              ||\n");
    printf("||                                                                              ||\n");
    printf("|+------------------------------------------------------------------------------+|\n");
    printf("+--------------------------------------------------------------------------------+\n");
    printf("\n\n\n");
  }

  //MPI_Barrier(MPI_COMM_WORLD);  printf("FINALIZE CALL----------pr=%d\n",myid);  MPI_Finalize(); return 0; //TERMINATION COMMANDS


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
  /*!
    \brief Initialize the 1D Vandermonde Matrix, V_{ij} = phi_j(r_i);
  */  
  V1D.resize(r.size1(),order+1);
  fullMatrix<scalar> P;
  for(int j=0;j<order+1;j++){
    JacobiP(r, 0, 0, j, P);
    for(int i=0;i<P.size1();i++) V1D(i,j) = P(i,0);
  }
}

void monovandermonde1d(const int order, const fullMatrix<double> r, fullMatrix<scalar> &V1D){
  /*!
    \brief Initialize the 1D Vandermonde Matrix, V_{ij} = (r_i)^j/factorial(j);
  */
  
  V1D.resize(r.size1(),order+1);
  for(int j=0;j<order+1;j++){
    for(int i=0;i<r.size1();i++){
      V1D(i,j) = pow(r(i,0),j)/(scalar)factorial(j);
    }
  }
}

void monovandermonde2d(const int order, const fullMatrix<double> r, fullMatrix<scalar> &V2D){
  /*!
    \brief Initialize the 2D Vandermonde Matrix, V = (x_i)^nx (y_i)^ny / (factorial(nx)*factorial(ny)).
    ith line = [1, x_i, y_i, x_i^2/2!, x_i*y_i, y_i^2/2!, x_i^3/3!, x_i^2*y_i/2!, x_i*y_i^2/2!, y_i^3/3!, ...]
  */

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

void LagMono2DTransforms(const int N_E, const int N_s, const int order, const int L2Msize1, const int L2Msize2, std::string ElemType, const fullMatrix<scalar> XYZNodes, const fullMatrix<scalar> XYZCen, fullMatrix<scalar> &Lag2Mono, fullMatrix<scalar> &Mono2Lag){
  /*!
    \brief Returns the transforms (and inverse) from Lagrange to Taylor polynomials.
    NB: In >< to the 1D transforms, these are in the physical space! So there is one transform per element.
  */
  
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


void getPowersXYZG(const int N_E, const int N_s, const int N_G, const int N_N, const int M_B, const int order, const fullMatrix<scalar> XYZG, const fullMatrix<scalar> XYZCen, const int* neighbors, const fullMatrix<scalar> shifts, scalar* powers){
  /*!
    \brief Get the powers of XZYG-XYZCen for each element and his neighbors/
    This is precalculated for increased speed in 2D limiting.
  */
    
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


int getTaylorDerIdx2DLength(const int order){
  /*! Calculate the length of these indexes*/
  int L = 0;
  for(int p = order; p >= 0; p--) L+=(p+1)*(p+2)/2;  // DOES THIS WORK ONLY FOR TRIANGLES?
  return L;
}

void getTaylorDerIdx2D(const int order, int* TaylorDxIdx, int* TaylorDyIdx){
  /*!
    \brief Get the x and y derivative index.
    Basically gives you the index for various derivatives of a Taylor
    polynomial. order is the DG order dxIdx =[idx for 0th derivative wrt
    x, idx for 1st derivative wrt x, ...]
  */

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

  /*!
    \brief Permutations for the cartesian mesh.
    This function should be used for a cartesian mesh.  The idea is to
    find the permutation matrix to go from the numbering system in gmsh
    to a system where the nodes are numbered in increasing order for
    increasing x and decreasing y. Best shown by example: for p=2, the
    values of U are stored at the following indexes:

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
  /*!
    \brief This function returns the 1D transforms for 2D transforms of Lagrangian basis.
    Ax = Lag2MonoX * U, Ay = MonoX2MonoY * Ax, U = MonoY2Lag Ay.
    This might be unnecessarily complicated (bc of the use of the GL
    points) or very smart (i.e stable bc of the use of GL points)
  */

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
