/*!
  \file mixedform.cc for dealing with 2nd order derivatives
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup gradient
*/

#include "mixedform.h"

void COPYCOPY_get_element_types(const int order, int &msh_hex, int &msh_qua, int &msh_tri, int &msh_lin){
  if      (order==0)  {msh_hex = MSH_HEX_8; msh_qua = MSH_QUA_4;    msh_tri = MSH_TRI_3;    msh_lin = MSH_LIN_2;  }
  else if (order==1)  {msh_hex = MSH_HEX_8; msh_qua = MSH_QUA_4;    msh_tri = MSH_TRI_3;    msh_lin = MSH_LIN_2;  }
  else if (order==2)  {msh_hex = MSH_HEX_27; msh_qua = MSH_QUA_9;    msh_tri = MSH_TRI_6;    msh_lin = MSH_LIN_3;  }
  else if (order==3)  {msh_hex = MSH_HEX_64; msh_qua = MSH_QUA_16;   msh_tri = MSH_TRI_10;   msh_lin = MSH_LIN_4;  }
  else if (order==4)  {msh_hex = MSH_HEX_125; msh_qua = MSH_QUA_25;   msh_tri = MSH_TRI_15;   msh_lin = MSH_LIN_5;  }
  else if (order==5)  {/*msh_hex = MSH_HEX_216;*/ msh_qua = MSH_QUA_36;   msh_tri = MSH_TRI_21;   msh_lin = MSH_LIN_6;  }
  else if (order==6)  {/*msh_hex = MSH_HEX_343; */msh_qua = MSH_QUA_49;   msh_tri = MSH_TRI_28;   msh_lin = MSH_LIN_7;  }
  else if (order==7)  {/*msh_hex = MSH_HEX_512; */msh_qua = MSH_QUA_64;   msh_tri = MSH_TRI_36;   msh_lin = MSH_LIN_8;  }
  else if (order==8)  {/*msh_hex = MSH_HEX_729; */msh_qua = MSH_QUA_81;   msh_tri = MSH_TRI_45;   msh_lin = MSH_LIN_9;  }
  else if (order==9)  {/*msh_hex = MSH_HEX_1000; */msh_qua = MSH_QUA_100;  msh_tri = MSH_TRI_55;   msh_lin = MSH_LIN_10; }
  else if (order==10) {/*msh_hex = MSH_HEX_1331; */msh_qua = MSH_QUA_121;  msh_tri = MSH_TRI_66;   msh_lin = MSH_LIN_11; }
  else {printf("Invalid order number.");}
}


void SurfForSigma(deck inputs, int N_s, int p, int om, int M_G, int N_N, int M_s, int a, /*int* rog_invmap_local,*/int* invmap, int* place, fullMatrix<double> pointsF, fullMatrix<double> weightF, scalar* JFace, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> OutwardNormal, fullMatrix<scalar> psi_w, /*int* indices_support, *//*scalar* Sa_return*/ fullMatrix<scalar> &Sa_return)
{
  /*!
    \brief populate an element's surface term to solve for the auxiliary variable.
    Actually, just one of the D submatrices of that matrix, whichever a is 
    specified in SurfForSigma			
    \param[in] inputs information from the deck
    \param[in] N_s solution points per element
    \param[in] p solution order
    \param[in] om = the element we are working on
    \param[in] M_G = quadrature points per interface
    \param[in] N_N neighbors per element
    \param[in] M_s supported nodes per face of reference element
    \param[in] a derivative direction
    \param[in] invmap gets the element index addresses for the supported basis functions along a face
    \param[in[ pointsF = (D-1) reference geometry: point for 1D, line for 2D, square/simplex for 3D
    \param[in] weightF = quadrature weights on the D-1 reference geometery
    \param[in] JFace = detJ along the face for face quadrature
    \param[in] XYZNodes the physical coordinates of solution nodes
    \param[in] OutwardNorrmal = outward normal from the element on each interface. Again, constant on each interface
    input: indices_support should have M_s=(p+1)^(D-1) entries, each of which is an element basis function index that is supported on the given interface
    \param[in] psi_w face-supported DG basis functions multiplied by face quadrature weights
    \param[out] Sa_return is a directional component of the surface residual matrix for mixed-form approach
*/

  //Outward normal and JaCFace: Take these from the parent routine
  
  //Build the basis functions along the cell interfaces (this is psi, but psi_w has already combined with quadr weights)
  int verbose = 0;
  
  //Assemble each interface's block in to the output matrix
  fullMatrix<scalar> output (N_s, N_N*M_G);
  if (verbose == 1){ printf("Entered SurfForSigma: om=%d, a=%d\n", om,a);} 
  for (int H = 0; H < N_N; H++) //H is the local face address
    {
      if (verbose == 1){  printf("Working local face %d\n", H);}

      //Now working on a specific interface of the local element.
      fullMatrix<scalar> subBlock (N_s, M_G);
      for (int row = 0; row  < N_s; row++)
	{
	  for (int g = 0; g < M_G; g++)
	    {
	      subBlock(row,g) = 0.0;
	    }
	}
      //only the face-supported basis functions get to contribute:
      for (int m = 0; m < M_s; m++)
	{
	  //int row = invmap[((om*N_F+0)*M_s*N_N +  place[H]*M_s + m)*2+0]; //solIDX of the intefacre-supported basis function
	  int row = invmap[((om*N_F+0)*M_s*N_N +  H*M_s + m)*2+0]; //solIDX of the intefacre-supported basis function
	  //int row = rog_invmap_local[H*M_s + m]; //solIDX of the intefacre-supported basis function
	  if (verbose == 1)
	    {
	      printf("H=%d, m=%d: basis index = %d, Outward normal[a]=%f, JFace=%f\n", H,m, row,OutwardNormal(H,a), JFace[H]);
#ifdef TWOD
	      printf("\t\tNode Location: (%f,%f)\n", XYZNodes(row, om*D + 0), XYZNodes(row, om*D + 1));
#endif
#ifdef THREED
	      printf("\t\tNode Location: (%f,%f,%f)\n", XYZNodes(row, om*D + 0), XYZNodes(row, om*D + 1), XYZNodes(row, om*D + 2));
#endif
	    }
	  for (int g = 0; g < M_G; g++)
	    {
	      //psi_w contains quadrature weight multiplied by face-supported basis function
	      //subBlock(row, g) += psi_w(g, m) * OutwardNormal(H,a) * JFace[H];
	      //subBlock(row, g) += psi_w(g, m) * OutwardNormal(H,a) * JFace[H];
	      //	      subBlock(row, g) += psi_w(g, m) * OutwardNormal(H,a) * JFace[H];
	      //09/20/2017: I don't know why, but the gradient system requires
	      //square of JFace in 3D to give proper behavior. Using power function
	      //to maintain good behavior in 1D and 2D cases

	      //subBlock(row, g) += psi_w(g, m) * OutwardNormal(H,a) * pow(JFace[H], D-1);
	      //okay, actually I decided to achieve congruence by messing with how JFace
	      //is calculated in the dg_functions/dg_jacobians_face routine.
	      subBlock(row, g) += psi_w(g, m) * OutwardNormal(H,a) * JFace[H];
	    }
	}
      if (verbose == 1){      printf("subBlock for interface %d: outward normal in a direction is %f\n", H,OutwardNormal(H,a));
      for (int row = 0; row < N_s; row++)
	{
	  printf("row %d: ",row);
	  for (int g = 0; g < M_G; g++)
	    {
	      printf("%f, ",subBlock(row,g));
	    }
	  printf("\n");
	}
      }
      /*
      for (int row = 0; row < N_s; row++)
	{
	  for (int g = 0; g < M_G; g++)
	    {
	      subBlock(row,g) = OutwardNormal(H,a) * weightF(g,0) * phi_trace(H*M_G + g , row) * JFace[H];
	    }
	} //end row loop for populating subBlock
      */
      //Send the local face's subBlock to the master output matrix
      for (int row = 0; row < N_s; row++)
	{
	  for (int g = 0; g < M_G; g++) //quadrature node on the interface
	    {
	      output(row, H*M_G + g) = subBlock(row, g);
	    }
	} //end the row loop for relaying subBlock to output

    } //end the gig H loop
  
  //Superfluous step, such that "output" is not in function declaration:
  //Sa_return is the matrix this routine was built to calculate
  Sa_return = output;
  /*
  //with output matrix populated, serialize it
  for (int row = 0; row < N_s; row++)
    {
      for (int col = 0; col < N_N*M_G; col++)
	{
	  int slot = row*N_N*M_G + col;
	  Sa_return[slot] = output(row,col);
	}
    }
  */
  //Sa_return is the value of interest, this subroutine is concluded
}

void VolForSigma(int om, deck inputs, int QuadRule, int N_G, int M_G, int N_s, int p, int a, scalar* invJac, scalar detJ, scalar* detJ_full, /*scalar* Va_return*/ fullMatrix<scalar> &Va_return)
{
  /*!
    \brief: Generate the volume residual matrix corresponding to a specficif directional component
    in the auxiliary weak form, perhaps the term 'mixed-form' is more familiar.
    \param[in] om local element address
    \param[in] inputs is the stuff from the deck, I want it so I can build basis functions
    \param[in] QuadRule quadrature precision demanded by main.cc
    \param[in] N_G quadrature points per element
    \param[in] M_G quadrature points per interface
    \param[in] N_s is nodes pere elemenet
    \param[in] p is element polynomial order
    \param[in] a direction in which derivative is being calculated
    \param[in] invJac is full inverse jacobian matrix, necessary for gradPhi
    \param[in] detJ is determinant of the element's jacobain, to perform integration
    \param[in] detJ_full is determinant of the element's jacobain, to perform integration
    \param[out] Va_return contribution of element's DOF to local sigma approximation in direction a
  */
  //Hint on matrix multiplication:
  /*
    Recovery matrix is LHS^-1 * RHS
    found the matrix multiplication command in main.cc, subroutine  LagMono2DTransformsCartesian, see last line
    fullMatrix<scalar> MatReco(N_basis, 2*N_s);
    MatReco.gemm(LHS_inv, RHS);
  */
  /*
    10/31/2017: I want to import detJ from main.cc, so this subroutine must use
    same quadrature precision as main.cc
   */
  int verbose = 0;
  if (verbose == 1){printf("Entering VolForSigma: om=%d, a=%d, element's detJ=%f\n",om,a,detJ);}
  fullMatrix<scalar> Beta(N_s, N_G);
  fullMatrix<scalar> Phi(N_G, N_s);
  //Phi is the easy matrix, do that first
  //Maybe one day, I'll want to use a higher resolution quadrature rule
  //for this procedure. So, I'm reforming the solution basis
  //here instead of importing it.

  //1) Get reference basis on reference element at all N_G quadrature points,
  //Quite a bit of this code is copied from early segments of main.cc because
  //I need to populate the DG basis functions
  // Get the method order
  int order = inputs.getOrder();
  bool order0 = false; if (order==0) {order0 = true; order = 1;}
  int elem_type;
  int face_type;
  int msh_hex;
  int msh_qua;
  int msh_tri;
  int msh_lin;
  int nsides; // this offsets j in buildInterfaces function
  int N_N;    // number of neighbors to an element
  scalar refArea; // area of reference element
  COPYCOPY_get_element_types(order, msh_hex, msh_qua, msh_tri, msh_lin);
  if     (inputs.getElemType() == "lin"){face_type = MSH_PNT, elem_type = msh_lin; nsides = 0; N_N = 2;}
  else if(inputs.getElemType() == "tri"){face_type = msh_lin, elem_type = msh_tri; nsides = 3; N_N = 3; refArea = 0.5;}
  else if(inputs.getElemType() == "qua"){face_type = msh_lin, elem_type = msh_qua; nsides = 4; N_N = 4; refArea = 4;} 
  else if(inputs.getElemType() == "hex"){face_type = msh_qua, elem_type = msh_hex; nsides = 6; N_N = 6; refArea = 8;}
  else printf("Invalid element type in deck");
  const polynomialBasis *basis  = polynomialBases::find(elem_type);  // for the element
  const std::vector<std::vector<int> > &closures = basis->closures;

  fullMatrix<double> points, weight;
  //FOR IRREGULAR ELEMENT SHAPES THIS QUADRATURE ORDER NEEDS TO BE HIGHER
  if     (inputs.getElemType() == "lin") gaussIntegration::getLine(QuadRule, points, weight);
  else if(inputs.getElemType() == "tri") gaussIntegration::getTriangle(QuadRule, points, weight);
  else if(inputs.getElemType() == "qua") gaussIntegration::getQuad(QuadRule, points, weight);
  else if(inputs.getElemType() == "hex") gaussIntegration::getHexahedron(QuadRule, points, weight);

  // fullMatrix<scalar> Phi (N_G , N_s); 
  //  printf("Building the basis matrix locally, N_G=%d, N_s=%d\n", N_G, N_s);
  fullMatrix<double> phiD (N_G , N_s); 
  basis->f (points, phiD);
  for(int g = 0; g < N_G; g++){
    for(int i = 0; i < N_s; i++){
      Phi(g,i) = (scalar)phiD(g,i);
      //    printf("Phi(%d,%d)=%f\n",g,i,Phi(g,i));
    }
  } 

  //Next: the beta integration matrix.
  //Now, this necessitates basis gradient wrt physical coordinate a, so we have a headache:
  /*
    //1D jacobian interpretation:
    dxi_dx = invJac[e*N_G*D*D + 0 + 0 + 0];

    //2D jacobian interpretation:
    scalar dxi_dx  = invJac[e*N_G*D*D + g*D*D + 0 + 0];
    scalar deta_dx = invJac[e*N_G*D*D + g*D*D + 0 + 1];
    scalar dxi_dy  = invJac[e*N_G*D*D + g*D*D + D + 0];
    scalar deta_dy = invJac[e*N_G*D*D + g*D*D + D + 1];
  */

  if (verbose == 1)
    {
      printf("The local elements of invJac\n");
      for (int g = 0; g < N_G; g++)
	{
	  printf("--quadrature point %d:---\n", g);
	  for (int a1 = 0; a1 < D; a1++){
	    printf("\trow %d: ",a1);
	    for (int a2 = 0; a2 < D; a2++){
	      printf("%f, ",invJac[om*N_G*D*D + g*D*D + a1*D + a2]);
	    }
	    printf(")\n");
	  }
	  printf("\n");
	}
    }
  // printf("Grabbing basis derivatives wrt reference coordinates\n");
  double grads[N_s][3];
  fullMatrix<scalar> dphi(N_G*D, N_s);
  for(int g = 0; g < N_G; g++){
    basis->df(points(g,0),points(g,1),points(g,2),grads);
    for(int alpha = 0; alpha < D; alpha ++){
      //this alpha is the reference coordinate, beware
      for(int i = 0; i < N_s; i++){
  	dphi(g*D+alpha,i) = (scalar)grads[i][alpha];  //see paper for indexing p.6
      }	  
    }  
  }
  if (verbose == 1){
  for (int k = 0; k < N_s; k++)
    {
      printf("---node index %d:---\n", k);
      for (int g = 0; g < N_G; g++)
	{
	  printf("g=%d, xi = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",points(g,a));
	    }
	  printf(")\n");
	  printf("gradient at this node wrt ref coordinates is (");
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",dphi(g*D + a, k));
	    }
	  printf(")\n");
	}
      printf("\n");
    }
  }

  //Now, get the gradient wrt X_a at every quadrature point:
  //a is an input integer
  fullMatrix<scalar> dphi_dXa (N_s, N_G);
  for (int k = 0; k < N_s; k++)
    {
      for (int g = 0; g < N_G; g++)
	{
	  dphi_dXa(k,g) = 0.0;
	  for (int alpha = 0; alpha < D; alpha++)
	    {
	      //Is this proper usage of invJac? who knows.
	      /*
		More analysis, 10/23/2017: For quite some time, I thought that it was
		correct to go (a*D + alpha), but today, I disovered that the code was 
		really bad at the auxiliary solve on simplex elements. So, I changed
		to (alpha*D+a) here and results look much better.
	       */
	      //dphi_dXa(k,g) += dphi(g*D+alpha, k) * invJac[om*N_G*D*D + g*D*D + a*D + alpha];
	      dphi_dXa(k,g) += dphi(g*D+alpha, k) * invJac[om*N_G*D*D + g*D*D + alpha*D + a];
	      //dphi_dXa(k,g) += dphi(g*D+alpha, k) * invJac[om*N_G*D*D + g*D*D + alpha*D + a];
	    }
	}
    }
  if (verbose == 1)
    {
      printf("dphi_dXa:\n");
      for (int k = 0; k < N_s; k++)
	{
	  printf("Basis index %d\n",k);
	  for (int g = 0; g < N_G; g++)
	    {
	      printf("derivative = %f\n",dphi_dXa(k,g));
	    }
	}
    
    }
  //I now have physical gradient at all N_G quadrature points. Also have quadrature weights.
  for (int row = 0; row < N_s; row++)
    {
      for (int g = 0; g < N_G; g++)
	{
	  //Beta(row, g) = dphi_dXa(row, g) * weight(g,0) * detJ;
	  //Use accurate detJ stturcture:
	  Beta(row, g) = dphi_dXa(row, g) * weight(g,0) * detJ_full[om*N_G + g];
	}
    }
  /*
  printf("Built the Beta integration matrix\n");
  for (int row = 0; row < N_s; row++)
    {
      printf("row %d:\t\t",row);
      for (int col = 0; col < N_G; col++)
	{
	  printf("%4.3f, ", Beta(row,col));
	}
      printf("\n");
    }
  */
  //multiply the volumne integration matrix by basis function matrix
  //to get the necessary volume matrix:
  fullMatrix<scalar> OUT(N_s, N_s);
  OUT.gemm(Beta, Phi);
  Va_return = OUT;
  //  printf("Multiplued Beta x Phi to get the output matrix\n");
  /*
  //scalar* Va_return = new scalar(N_s, N_s);
  //with output matrix populated, serialize it
  for (int row = 0; row < N_s; row++)
    {
      for (int col = 0; col < N_s; col++)
	{
	  int slot = row*N_s + col;
	  Va_return[slot] = OUT(row,col);
	}
    }
  */
  //Va_return is necessary output, subroutine concluded
} 

void BuildSigMatrices(deck inputs, int QuadRule, int N_E, int N_s, int p, int N_G, int M_G, int N_N, int M_s, int* invmap, fullMatrix<double> &pointsF, fullMatrix<double> &weightF, fullMatrix<scalar> &JF, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> normals, scalar* h_invJac, fullMatrix<scalar> &J, int* Alt_FaceFromElem, int* BR2_map, fullMatrix<scalar> phi, scalar* Minv_Master, fullMatrix<scalar> psi_w, scalar* detJ_full, scalar* serial_SigSurf, scalar* serial_SigVol)
{
  /*
  for (int t = 0; t < JF.size1(); t++)
    {
      printf("JF(%d) = %f\n",t,JF(t,0));  
    }
  */
  /*!
    \brief: Subroutine to build the matrices that get auxiliary variable from UCommon and Uhat
    I'm accounting for surface and volume contributions separately to 
    eventually lower communication cost in the code.
    \param[in] inputs is the master code parameters from the deck
    \param[in] QuadRule quadrature precision from main.cc
    \param[in] N_s is DOF per element
    \param[in] p = polynomial order
    \param[in] N_G = quadrature points per element interior
    \param[in] M_G = quadrature points per interface, one-sided
    \param[in] M_s = supported basis functions per face of one element (not counting the supported functions from the neighbor element)
    \param[in] invmap: structure that yields the basis function indices corresponding to the M_s suppoerted functions on a given face of an element
    \param[in] pointsF: (D-1)-dimensional set of quadrature points on reference interface
    \param[in] weightF: quadrature weights for reference interface
    \param[in] JF = jacobians of the interfaces
    \param[in] XYZNodes = all the solution node locations, all elements
    \param[in] normals = interface normals. I think these are outward from A-element, not sure
    \param[in] invJac = inverse Jacobian matrix (not just determinant) for each element.
    \param[in] J = det|Jacobian| for each element
    \param[in] detJ_full a reliable det(Jacobain) for each element
    
    \param[out] serial_SigSurf gets the contribution of interface common U values to gradient approximation
    \param[out[ serial_SigVol gets the contribution of element DOF to gradient approximation.
   */
  int verbose = 0;
  //Before anything else: reorganize invmap.
  fullMatrix<scalar> rog_invmap(N_E, M_s*N_N);
  /*
  int order_t[N_E][N_N];
  for (int om = 0; om < N_E; om++)
    {
      //      int order_t[N_N];
      //start: assume that each t in FaceIndex is the lowest t
      for (int H = 0; H < N_N; H++)
	{
	  order_t[om][H] = 0;
	}
      for (int H = 0; H < N_N; H++)
	{
	  int t_local = FaceIndex[H];
	}
      for (int H2 = 0; H2 < N_N; H2++)
	{
	  if (H != H2 && t_local > FaceIndex[H])
	    {
	      order_t[om][H] +=1;
	    }
	}
    }

  for (int om = 0; om < 1; om++)
    {
      for (int H = 0; H < N_N; H++)
	{
	  printf("H=%d, tGlo = %d, place = %d\n", H, FaceIndex[H], order_t[om][H]);
	}
    }
  */
  int yes = 1;
  if (yes == 1)
    {
      for (int om = 0; om < N_E; om++) //Run through all elements, set up the auxiliary solve
	{
	  if (verbose == 1){     printf("In BuildSigMatrices, om=%d\n",om);}
	  //Output matrices:
	  fullMatrix<scalar> FullSurf (D*N_G, N_N*M_G);
	  fullMatrix<scalar> FullVol (D*N_G, N_s);
	  //these output matrices are stackeed; there are D substructures
	  //in each, and each super-row corresponds to physical gradient in different direction
	  fullMatrix<scalar> Minv_local (N_s, N_s);
	  //Extract the element's inverse mass matrix from global storage
	  //copied this form straight from where Marc builds Minv in dg_functions.h
	  for(int k = 0; k < N_s; k++){
	    for(int i = 0; i < N_s; i++){
	      //	Minv[(e*N_s+k)*N_s+i] = M_inv(i,k);
	      //Sinnce Minv_master is column-major sorted, must
	      //use non-intuitive indexing to extract the local SSAM.
	      //Validated by eperimentation: SSAM should be its own transpose.
	      //Actually, doesn't matter for SSAM, will need to experiment elsewhere
	      Minv_local(i,k) = Minv_Master[(om*N_s+k)*N_s+i];
	      //Minv_local(i,k) = Minv_Master[(om*N_s+i)*N_s+k];
	    }
	  }
	  if (verbose == 1){
	    printf("Local inverse mass matrix:\n");
	    for (int row = 0; row < N_s; row++)
	    {
	    printf("row %d: ",row);
	    for (int col = 0; col < N_s; col++)
	    {
	    printf("%f, ",Minv_local(row,col));
	    }
	    printf("\n");
	    }
	    }
	  //Gonna need a bit of information about this element :)
	  //JFace is just the face jacobians of the local element. Must extract from JF
	  scalar* JFace = new scalar[N_N];
	  
	  //face jacobians are stored in face indexing, so I need to know the global interface
	  //indices corresponding to the element
	  int FaceIndex[N_N];
	  for (int j = 0; j < N_N; j++)
	    {
	      //FaceIndex[j] = Alt_FaceFromElem[om*N_N + j];
	      FaceIndex[j] = Alt_FaceFromElem[om*(N_N) + j];
	    }
	  for (int H = 0; H < N_N; H++)
	    {
	      //  printf("JF[%d] = %f\n", FaceIndex[H]*2 + 0, JF(FaceIndex[H]*2 + 0, 0));
	      JFace[H] = JF(FaceIndex[H]*2 + 0, 0); //JF has two rows per interface, not sure why.
	      //Actually, I think it's two values because there's one for each interface.
	      //They are identical, but it's good to follow the form of duplicate
	      //storage for both sides of interface.
	    }
	  if (verbose == 1)
	    {
	      printf("FaceIndex and corresponding face jacobians:\n");
	      for (int H = 0; H < N_N; H++)
		{
		  printf("H=%d, tGlo = %d, JFace=%f, omA=%d, omB=%d\n", H, FaceIndex[H], JFace[H], BR2_map[FaceIndex[H]*4 + 2*0 + 0], BR2_map[FaceIndex[H]*4 + 2*1 + 0]);
		}
	    }
	  int order_t[N_N];
	  //start: assume that each t in FaceIndex is the lowest t
	  for (int H = 0; H < N_N; H++)
	    {
	      order_t[H] = 0;
	    }
	  for (int H = 0; H < N_N; H++)
	    {
	      int t_local = FaceIndex[H];
	
	      for (int H2 = 0; H2 < N_N; H2++)
		{
		  if (H != H2 && t_local > FaceIndex[H2])
		    {
		      order_t[H] +=1;
		    }
		}
	    }

	  /*
	    Need to talk about why we have order_t. In the SurfForSigma
	    subroutine, I use Marc's invmap structure. Now, for a given
	    element, invmap's entries give the element's DOF indices
	    corresponding to each support point on each interface.
	    However, invmap is not ordered according to the local
	    face indexing. Instead, the lowest tGlo value goes first,
	    and the highet tGlo value goes last. Now, since I'm 
	    using the BR2_map array to get the interface
	    normals, I need invmap to give information about
	    node mappings in the proper order. So, when
	    looking for supported DOF indices on an interface
	    in SurfForSigma, instead of using H, I call place[H]
	    to grab the entries of invmap corresponding to
	    local face H.

	    05/23/2017: Trying to get around this issue by instead
	    creating Alt_FaceFromElem, which stores an element's
	    global interface addresses in the order
	    that they appear in invmap
	  */



	  //With the order of the tGlo addresses known, populate the
	  //local element's reorganized inverse map.
	  int rog_invmap_local[N_N*M_s];
	  for (int H = 0; H < N_N; H++)
	    {
	      for (int m = 0; m < M_s; m++)
		{
		  rog_invmap_local[H*M_s + m] = invmap[((om*N_F+0)*M_s*N_N +  order_t[H]*M_s + m)*2+0]; //solIDX of the intefacre-supported basis function
		  rog_invmap(om, H*M_s + m) = rog_invmap_local[H*M_s + m];
		}
	    }


	  scalar detJ_local = J(om,0);
	  //The normal is tricky. Since it changes signs depending on the outward normal,
	  //I must be sure that the one I'm picking up corresponds to the local element.
	  //Can take care of this thru BR2_map
	  fullMatrix<scalar> OutwardNormal (N_N,D);
	  for (int H = 0; H < N_N; H++)
	    {
	      int tFace = FaceIndex[H];
	      //	  printf("H=%d, tFace=%d\n", H, tFace);
	      //I think if sig_index ib BR2_Map is as big or larger than N_N, normal is reversed
	      int omA_face = BR2_map[tFace*4 + 2*0 + 0];
	      int omB_face = BR2_map[tFace*4 + 2*1 + 0];
	      scalar mult_normal = 1.0;
	      if (om == omA_face)
		{
		  //do nothing, normal is pointing in proper direction
		}
	      else if (om == omB_face)
		{
		  mult_normal = -1.0; //reverse normal, because we are in B element
		}
	      else
		{
		  printf("CATASTROPHE in mixedform.cc, BuildSigMatrices!!! BR2_Map has a mismatch\n");
		}
	      for (int a = 0; a < D; a++)
		{
		  //h_normals is in column-major form, need to massage it a bit:
		  //	      OutwardNormal(H,a) = h_normals[tFace*D + a]*mult_normal;
		  //OutwardNormal(H,a) = h_normals[a*M_T + tFace]*mult_normal;
		  //For the organization of fullMatrix normals, see buildNormals in simplemesh.cc
		  OutwardNormal(H,a) = normals(a, tFace) * mult_normal;
		  //OutwardNormal(order_t[H],a) = normals(a, tFace);
		}
	    }


	  if (verbose == 1)
	    {
	    printf("Got the outwardNormal on each face:\n");
	    for (int H = 0; H < N_N; H++)
	    {
	    printf("Face %d, tGlo=%d  ",H, FaceIndex[H]);
	    for (int a = 0; a < D; a++)
	    {
	    printf("n_(%d) = %f,  ",a, OutwardNormal(H,a));
	    }
	    printf("\n");
	    }
	    }
	  for (int a = 0; a < D; a++) //directional component of sigma that we need
	    {
	      fullMatrix<scalar> DirSurf_weak (N_s, N_N*M_G);
	      fullMatrix<scalar> DirVol_weak (N_s, N_s);
	      //scalar* unrolled_DirSurf_weak = new scalar[N_s*N_N*M_G];
	      //scalar* unrolled_DirVol_weak  = new scalar[N_s*N_s];
	      //scalar* random_matrix;
	      // delete[] unrolled_DirSurf_weak;
	      //delete[] unroll
	      /*
		DirSurf_weak and DirVol_weak are the matrices, that form
		the auxiliary weak form. To solve for SigHat, would have to multiply
		by inverse mass matrix (we'll get there)
	      */
	  
	      if (verbose == 1){  printf("Direction = %d: Calling SurfForSigma\n",a);}

	      SurfForSigma(inputs, N_s, p, om, M_G, N_N, M_s, a, invmap, order_t, pointsF, weightF, JFace, XYZNodes, OutwardNormal, psi_w, /*unrolled_DirSurf_weak*/ DirSurf_weak);
	      //SurfForSigma(inputs, N_s, p, om, M_G, N_N, M_s, a, invmap, pointsF, weightF, JFace, XYZNodes, OutwardNormal, psi_w, /*unrolled_DirSurf_weak*/ DirSurf_weak);
	      //SurfForSigma(inputs, N_s, p, om, M_G, N_N, M_s, a, rog_invmap_local, pointsF, weightF, JFace, XYZNodes, OutwardNormal, psi_w, /*unrolled_DirSurf_weak*/ DirSurf_weak);

	      //  printf("Back in SigSurfBuild\n");
	      // printf("Direction = %d: Calling VolForSigma\n",a);
	      VolForSigma(om, inputs, QuadRule, N_G, M_G, N_s, p, a, h_invJac, detJ_local/*,random_matrix, unrolled_DirVol_weak*/, detJ_full, DirVol_weak);
	  
	      //  printf("Back in SigSurfBuild\n");

	      /*
	      //Need to put these unrolled matrices back in matrix form
	      for (int row = 0; row < N_s; row++)
	      {
	      for (int col = 0; col < N_N*M_G; col++)
	      {
	      int slot = row*N_N*M_G + col;
	      DirSurf_weak(row,col) = unrolled_DirSurf_weak[slot];
	      }
	      }
	      for (int row = 0; row < N_s; row++)
	      {
	      for (int col = 0; col < N_s; col++)
	      {
	      int slot = row*N_s + N_s;
	      DirVol_weak(row,col) = unrolled_DirVol_weak[slot];
	      }
	      }
	      */

	      //Now, multiply each of these weak form matrices by inverse mass matrix
	      fullMatrix<scalar> SigHat_cont_Surf (N_s, N_N*M_G);
	      fullMatrix<scalar> SigHat_cont_Vol (N_s, N_s);
	      SigHat_cont_Surf.gemm(Minv_local, DirSurf_weak);
	      SigHat_cont_Vol.gemm(Minv_local, DirVol_weak);

	      //	  printf("om=%d, a=%d: performed SSAM multiplication\n",om,a);
	      if (verbose == 1)
		{
		  printf("Output from the SurfForSigma routine: DirSurfWeak\n");
		  for (int row = 0; row < N_s; row++)
		    {
		      printf("resi for row %d:\n", row);
		      for (int H = 0; H < N_N; H++)
			{
			  printf("|");
			  for (int g = 0; g < M_G; g++)
			    {
			      printf("%f,",DirSurf_weak(row, H*M_G+g));
			    }
			}
		      printf("\n");
		    }
		  printf("Output from the VolForSigma routine: DirVolWeak\n");
		  for (int row = 0; row < N_s; row++)
		    {
		      printf("resi for row %d:\n", row);
		      for (int col = 0; col < N_s; col++)
			{
			  printf("%f,",DirVol_weak(row, col));
			}
		      
		      printf("\n");
		    }
		  printf("After SSAM multiplication: SigHat_cont_Surf to get surface contribution to auxiliary coefficients:\n");
		  for (int row = 0; row < N_s; row++)
		    {
		      printf("coefficients for row %d:\n", row);
		      for (int H = 0; H < N_N; H++)
			{
			  printf("|");
			  for (int g = 0; g < M_G; g++)
			    {
			      printf("%f,",SigHat_cont_Surf(row, H*M_G+g));
			    }
			}
		      printf("\n");
		    }
		  printf("Sum of each row of SigHat_cont_Surf matrix\n");
		  for (int row = 0; row < N_s; row++)
		    {
		      scalar sumRow = 0;
		      for (int j = 0; j < N_N*M_G; j++)
			{
			  sumRow += SigHat_cont_Surf(row, j);
			}
		      printf("row %d: sum = %f\n",row,sumRow);
		    }
		  printf("Sum of each row of SigHat_cont_Vol matrix\n");
		  for (int row = 0; row < N_s; row++)
		    {
		      scalar sumRow = 0;
		      for (int j = 0; j < N_G; j++)
			{
			  sumRow += SigHat_cont_Vol(row, j);
			}
		      printf("row %d: sum = %f\n",row,sumRow);
		    }

		  
		  printf("The phi matrix:\n");
		  for (int row = 0; row < N_G; row++)
		    {
		      printf("g=%d: ",row);
		      for (int col = 0; col < N_s; col++)
			{
			  printf("%f, ",phi(row,col));
			}
		      printf("\n");
		    }
		}

	      /*
		Now, we have the matrices that can get the auxiliary coefficients
		for an element, in a given derivative direction, 
		given the Ucommon distribution and the elemen't Uhat collection.
		But that's not enough for me. Multiply by the element's shape functions
		to go directly to the gradient approximation at the element's quadrature points
	      */
	      fullMatrix<scalar> SigAct_cont_Surf (N_G, N_N*M_G);
	      fullMatrix<scalar> SigAct_cont_Vol (N_G, N_s);
	      SigAct_cont_Surf.gemm(phi, SigHat_cont_Surf);
	      SigAct_cont_Vol.gemm(phi, SigHat_cont_Vol);

	  
	  
	      //Okay, one more thing. The volume contribution needs a negative sign
	      for (int row = 0; row < N_G; row++)
		{
		  for (int col = 0; col < N_s; col++)
		    {
		      SigAct_cont_Vol(row,col) = -1.0 * SigAct_cont_Vol(row,col);
		    }
		}
	      /*
		printf("om=%d, a=%d: multiplied by phi to get the U->GradU transformation\n",om,a);
		printf("The derivative contribution of UhCommon:\n");
		for (int g = 0; g < N_G; g++)
		{
		printf("g=%d:\t",g);
		for (int i = 0; i < N_N*M_G; i++)
		{
		printf("%6.4f, ", SigAct_cont_Surf(g,i));
		}
		printf("\n");
		}
		printf("\n\nThe derivative contribution of Uhat\n");
		for (int g = 0; g < N_G; g++)
		{
		printf("g=%d\t",g);
		for (int i = 0; i < N_s; i++)
		{
		printf("%6.4f, ",SigAct_cont_Vol(g, i));
		}
		printf("\n");
		}
	      */
	      //05/22/2017: results up to here look belieavable for element 0. Mus test for more elem.

	      //A quick check: Plug in 1 for Uhcommon and Uhat everywhere, gradient should be zero
	      fullMatrix<scalar> Quiescent_Surf (N_N*M_G,1);
	      fullMatrix<scalar> Quiescent_Vol   (N_s, 1);
	      for (int row = 0; row < N_N*M_G; row++)
		{
		  Quiescent_Surf(row,0) = 1.0;
		}
	      for (int row = 0; row < N_s; row++)
		{
		  Quiescent_Vol(row,0) = 1.0;
		}
	      fullMatrix<scalar> VolCheck (N_G, 1);
	      fullMatrix<scalar> SurfCheck (N_G, 1);
	      VolCheck.gemm(SigAct_cont_Vol, Quiescent_Vol);
	      SurfCheck.gemm(SigAct_cont_Surf, Quiescent_Surf);
	      if (om == 0)
		{
		  //HOUSTON, WE HAVE A PROBLEm: THIS TEST IS NOT RETURNING ZERO-GRADIENT
		  fullMatrix<scalar> Check (N_G, 1);
		  //		  printf("The Checking matrix for stationary IC\n");
		  for (int row = 0; row < N_G; row++)
		    {
		      Check(row,0) = SurfCheck(row,0) + VolCheck(row,0);
		      //      printf("Check(%d,0) = %f\n",row,Check(row,0));
		    }
		}
	  

	      //Send these directional contributions to the output matrices
	      for (int row = 0; row < N_G; row++)
		{
		  for (int col = 0; col < N_N*M_G; col++)
		    {
		      FullSurf(a*N_G + row, col) = SigAct_cont_Surf(row, col);
		    }
		  for (int col = 0; col < N_s; col++)
		    {
		      FullVol(a*N_G + row, col) = SigAct_cont_Vol(row,col);
		    }
		}
	    } //end of direction 'a' loop

	  //Place FullSurf and FullVol in giant array format.
	  //I'd like to cast these directly in column-major form
	  //for a blasGemm call, but that may not be possible
	  //here because serial_SigSurf must be condensed
	  //down to M_T*M_G total columns
	  for (int row = 0; row < D*N_G; row++)
	    {
	      for (int col = 0; col < N_N*M_G; col++)
		{
		  serial_SigSurf[om*(D*N_G*N_N*M_G) + row*(N_N*M_G) + col] = FullSurf(row,col);
		}
	      for (int col = 0; col < N_s; col++)
		{
		  serial_SigVol[om*(D*N_G*N_s) + row*N_s + col] = FullVol(row,col);
		}
	    }
	  delete[] JFace;
	} //end of giant element loop
    }
  //serial_SigSurf and serial_SigVol can now be applied to approximate the gradient
  //over all element interiors. Subroutine concluded
  
}

void BuildAuxHatMatrices(deck inputs, int QuadRule, int N_E, int N_s, int p, int N_G, int M_G, int N_N, int M_s, int* invmap, fullMatrix<double> &pointsF, fullMatrix<double> &weightF, fullMatrix<scalar> &JF, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> normals, scalar* h_invJac, fullMatrix<scalar> &J, int* Alt_FaceFromElem, int* BR2_map, fullMatrix<scalar> phi, scalar* Minv_Master, fullMatrix<scalar> psi_w, scalar* detJ_full, scalar* serial_AuxHatSurf, scalar* serial_AuxHatVol)
{
  /*
  for (int t = 0; t < JF.size1(); t++)
    {
      printf("JF(%d) = %f\n",t,JF(t,0));  
    }
  */
  /*!
    \brief: Subroutine to build the matrices that get auxiliary variable from UCommon and Uhat
    I'm accounting for surface and volume contributions separately to 
    eventually lower communication cost in the code.
    \param[in] inputs is the master code parameters from the deck
    \param[in] QuadRule quadrature precision from main.cc
    \param[in] N_s is DOF per element
    \param[in] p = polynomial order
    \param[in] N_G = quadrature points per element interior
    \param[in] M_G = quadrature points per interface, one-sided
    \param[in] M_s = supported basis functions per face of one element (not counting the supported functions from the neighbor element)
    \param[in] invmap: structure that yields the basis function indices corresponding to the M_s suppoerted functions on a given face of an element
    \param[in] pointsF: (D-1)-dimensional set of quadrature points on reference interface
    \param[in] weightF: quadrature weights for reference interface
    \param[in] JF = jacobians of the interfaces
    \param[in] XYZNodes = all the solution node locations, all elements
    \param[in] normals = interface normals. I think these are outward from A-element, not sure
    \param[in] invJac = inverse Jacobian matrix (not just determinant) for each element.
    \param[in] J = det|Jacobian| for each element
    \param[in] detJ_full a reliable det(Jacobain) for each element
    
    \param[out] serial_AuxHatSurf gets the contribution of interface common U values to gradient approximation coefficients
    \param[out[ serial_AuxHatVol gets the contribution of element DOF to gradient approximation coefficients.
   */
  int verbose = 0;
  //Before anything else: reorganize invmap.
  fullMatrix<scalar> rog_invmap(N_E, M_s*N_N);
  /*
  int order_t[N_E][N_N];
  for (int om = 0; om < N_E; om++)
    {
      //      int order_t[N_N];
      //start: assume that each t in FaceIndex is the lowest t
      for (int H = 0; H < N_N; H++)
	{
	  order_t[om][H] = 0;
	}
      for (int H = 0; H < N_N; H++)
	{
	  int t_local = FaceIndex[H];
	}
      for (int H2 = 0; H2 < N_N; H2++)
	{
	  if (H != H2 && t_local > FaceIndex[H])
	    {
	      order_t[om][H] +=1;
	    }
	}
    }

  for (int om = 0; om < 1; om++)
    {
      for (int H = 0; H < N_N; H++)
	{
	  printf("H=%d, tGlo = %d, place = %d\n", H, FaceIndex[H], order_t[om][H]);
	}
    }
  */
  int yes = 1;
  if (yes == 1)
    {
      for (int om = 0; om < N_E; om++) //Run through all elements, set up the auxiliary solve
	{
	  if (verbose == 1){     printf("In BuildSigMatrices, om=%d\n",om);}
	  //Output matrices:
	  fullMatrix<scalar> FullSurf (D*N_G, N_N*M_G);
	  fullMatrix<scalar> FullVol (D*N_G, N_s);
	  //these output matrices are stackeed; there are D substructures
	  //in each, and each super-row corresponds to physical gradient in different direction
	  fullMatrix<scalar> Minv_local (N_s, N_s);
	  //Extract the element's inverse mass matrix from global storage
	  //copied this form straight from where Marc builds Minv in dg_functions.h
	  for(int k = 0; k < N_s; k++){
	    for(int i = 0; i < N_s; i++){
	      //	Minv[(e*N_s+k)*N_s+i] = M_inv(i,k);
	      //Sinnce Minv_master is column-major sorted, must
	      //use non-intuitive indexing to extract the local SSAM.
	      //Validated by eperimentation: SSAM should be its own transpose.
	      //Actually, doesn't matter for SSAM, will need to experiment elsewhere
	      Minv_local(i,k) = Minv_Master[(om*N_s+k)*N_s+i];
	      //Minv_local(i,k) = Minv_Master[(om*N_s+i)*N_s+k];
	    }
	  }
	  if (verbose == 1){
	    printf("Local inverse mass matrix:\n");
	    for (int row = 0; row < N_s; row++)
	    {
	    printf("row %d: ",row);
	    for (int col = 0; col < N_s; col++)
	    {
	    printf("%f, ",Minv_local(row,col));
	    }
	    printf("\n");
	    }
	    }
	  //Gonna need a bit of information about this element :)
	  //JFace is just the face jacobians of the local element. Must extract from JF
	  scalar* JFace = new scalar[N_N];
	  
	  //face jacobians are stored in face indexing, so I need to know the global interface
	  //indices corresponding to the element
	  int FaceIndex[N_N];
	  for (int j = 0; j < N_N; j++)
	    {
	      //FaceIndex[j] = FaceFromElem[om*N_N + j];
	      FaceIndex[j] = Alt_FaceFromElem[om*N_N + j];
	    }
	  for (int H = 0; H < N_N; H++)
	    {
	      //  printf("JF[%d] = %f\n", FaceIndex[H]*2 + 0, JF(FaceIndex[H]*2 + 0, 0));
	      JFace[H] = JF(FaceIndex[H]*2 + 0, 0); //JF has two rows per interface, not sure why.
	      //Actually, I think it's two values because there's one for each interface.
	      //They are identical, but it's good to follow the form of duplicate
	      //storage for both sides of interface.
	    }
	  if (verbose == 1)
	    {
	      
	      printf("FaceIndex and corresponding face jacobians:\n");
	      for (int H = 0; H < N_N; H++)
		{
		  printf("H=%d, tGlo = %d, JFace=%f\n", H, FaceIndex[H], JFace[H]);
		}
	    }
	  int order_t[N_N];
	  //start: assume that each t in FaceIndex is the lowest t
	  for (int H = 0; H < N_N; H++)
	    {
	      order_t[H] = 0;
	    }
	  for (int H = 0; H < N_N; H++)
	    {
	      int t_local = FaceIndex[H];
	
	      for (int H2 = 0; H2 < N_N; H2++)
		{
		  if (H != H2 && t_local > FaceIndex[H2])
		    {
		      order_t[H] +=1;
		    }
		}
	    }
	  /*
	    for (int H = 0; H < N_N; H++)
	    {
	    printf("H=%d, tGlo = %d, place = %d\n", H, FaceIndex[H], order_t[H]);
	    }

	    printf("The element's N_N*M_s invmap entries, straight from invmap with no doctoring:\n");
	    for (int i = 0; i < N_N*M_s; i++)
	    {
	    printf("entry %d: element node = %d\n",i, invmap[((om*N_F+0)*M_s*N_N + i)*2+0]);  
	    }
	  */

	  /*
	    Need to talk about why we have order_t. In the SurfForSigma
	    subroutine, I use Marc's invmap structure. Now, for a given
	    element, invmap's entries give the element's DOF indices
	    corresponding to each support point on each interface.
	    However, invmap is not ordered according to the local
	    face indexing. Instead, the lowest tGlo value goes first,
	    and the highet tGlo value goes last. Now, since I'm 
	    using the BR2_map array to get the interface
	    normals, I need invmap to give information about
	    node mappings in the proper order. So, when
	    looking for supported DOF indices on an interface
	    in SurfForSigma, instead of using H, I call place[H]
	    to grab the entries of invmap corresponding to
	    local face H.

	    05/23/2017: Trying to get around this issue by instead
	    creating Alt_FaceFromElem, which stores an element's
	    global interface addresses in the order
	    that they appear in invmap
	  */



	  //With the order of the tGlo addresses known, populate the
	  //local element's reorganized inverse map.
	  int rog_invmap_local[N_N*M_s];
	  for (int H = 0; H < N_N; H++)
	    {
	      for (int m = 0; m < M_s; m++)
		{
		  rog_invmap_local[H*M_s + m] = invmap[((om*N_F+0)*M_s*N_N +  order_t[H]*M_s + m)*2+0]; //solIDX of the intefacre-supported basis function
		  rog_invmap(om, H*M_s + m) = rog_invmap_local[H*M_s + m];
		}
	    }


	  scalar detJ_local = J(om,0);
	  //The normal is tricky. Since it changes signs depending on the outward normal,
	  //I must be sure that the one I'm picking up corresponds to the local element.
	  //Can take care of this thru BR2_map
	  fullMatrix<scalar> OutwardNormal (N_N,D);
	  for (int H = 0; H < N_N; H++)
	    {
	      int tFace = FaceIndex[H];
	      //	  printf("H=%d, tFace=%d\n", H, tFace);
	      //I think if sig_index ib BR2_Map is as big or larger than N_N, normal is reversed
	      int omA_face = BR2_map[tFace*4 + 2*0 + 0];
	      int omB_face = BR2_map[tFace*4 + 2*1 + 0];
	      scalar mult_normal = 1.0;
	      if (om == omA_face)
		{
		  //do nothing, normal is pointing in proper direction
		}
	      else if (om == omB_face)
		{
		  mult_normal = -1.0; //reverse normal, because we are in B element
		}
	      else
		{
		  printf("CATASTROPHE in mixedform.cc, BuildAuxHatMatrices!!! BR2_Map has a mismatch\n");
		}
	      for (int a = 0; a < D; a++)
		{
		  //h_normals is in column-major form, need to massage it a bit:
		  //	      OutwardNormal(H,a) = h_normals[tFace*D + a]*mult_normal;
		  //OutwardNormal(H,a) = h_normals[a*M_T + tFace]*mult_normal;
		  //For the organization of fullMatrix normals, see buildNormals in simplemesh.cc
		  OutwardNormal(H,a) = normals(a, tFace) * mult_normal;
		  //OutwardNormal(order_t[H],a) = normals(a, tFace);
		}
	    }


	  if (verbose == 1)
	    {
	    printf("Got the outwardNormal on each face:\n");
	    for (int H = 0; H < N_N; H++)
	    {
	    printf("Face %d, tGlo=%d  ",H, FaceIndex[H]);
	    for (int a = 0; a < D; a++)
	    {
	    printf("n_(%d) = %f,  ",a, OutwardNormal(H,a));
	    }
	    printf("\n");
	    }
	    }
	  for (int a = 0; a < D; a++) //directional component of sigma that we need
	    {
	      fullMatrix<scalar> DirSurf_weak (N_s, N_N*M_G);
	      fullMatrix<scalar> DirVol_weak (N_s, N_s);
	      //scalar* unrolled_DirSurf_weak = new scalar[N_s*N_N*M_G];
	      //scalar* unrolled_DirVol_weak  = new scalar[N_s*N_s];
	      //scalar* random_matrix;
	      // delete[] unrolled_DirSurf_weak;
	      //delete[] unroll
	      /*
		DirSurf_weak and DirVol_weak are the matrices, that form
		the auxiliary weak form. To solve for SigHat, would have to multiply
		by inverse mass matrix (we'll get there)
	      */
	  
	      if (verbose == 1){  printf("Direction = %d: Calling SurfForSigma\n",a);}

	      SurfForSigma(inputs, N_s, p, om, M_G, N_N, M_s, a, invmap, order_t, pointsF, weightF, JFace, XYZNodes, OutwardNormal, psi_w, /*unrolled_DirSurf_weak*/ DirSurf_weak);
	      //SurfForSigma(inputs, N_s, p, om, M_G, N_N, M_s, a, invmap, pointsF, weightF, JFace, XYZNodes, OutwardNormal, psi_w, /*unrolled_DirSurf_weak*/ DirSurf_weak);
	      //SurfForSigma(inputs, N_s, p, om, M_G, N_N, M_s, a, rog_invmap_local, pointsF, weightF, JFace, XYZNodes, OutwardNormal, psi_w, /*unrolled_DirSurf_weak*/ DirSurf_weak);

	      //  printf("Back in SigSurfBuild\n");
	      // printf("Direction = %d: Calling VolForSigma\n",a);
	      VolForSigma(om, inputs, QuadRule, N_G, M_G, N_s, p, a, h_invJac, detJ_local/*,random_matrix, unrolled_DirVol_weak*/, detJ_full, DirVol_weak);
	  
	      //  printf("Back in SigSurfBuild\n");

	      /*
	      //Need to put these unrolled matrices back in matrix form
	      for (int row = 0; row < N_s; row++)
	      {
	      for (int col = 0; col < N_N*M_G; col++)
	      {
	      int slot = row*N_N*M_G + col;
	      DirSurf_weak(row,col) = unrolled_DirSurf_weak[slot];
	      }
	      }
	      for (int row = 0; row < N_s; row++)
	      {
	      for (int col = 0; col < N_s; col++)
	      {
	      int slot = row*N_s + N_s;
	      DirVol_weak(row,col) = unrolled_DirVol_weak[slot];
	      }
	      }
	      */

	      //Now, multiply each of these weak form matrices by inverse mass matrix
	      fullMatrix<scalar> SigHat_cont_Surf (N_s, N_N*M_G);
	      fullMatrix<scalar> SigHat_cont_Vol (N_s, N_s);
	      SigHat_cont_Surf.gemm(Minv_local, DirSurf_weak);
	      SigHat_cont_Vol.gemm(Minv_local, DirVol_weak);


	      /*
		Now, we have the matrices that can get the auxiliary coefficients
		for an element, in a given derivative direction, 
		given the Ucommon distribution and the elemen't Uhat collection.
	      */
	  

	      //Send these directional contributions to the output matrices
	      /*
	      for (int row = 0; row < N_G; row++)
		{
		  for (int col = 0; col < N_N*M_G; col++)
		    {
		      FullSurf(a*N_G + row, col) = SigAct_cont_Surf(row, col);
		    }
		  for (int col = 0; col < N_s; col++)
		    {
		      FullVol(a*N_G + row, col) = SigAct_cont_Vol(row,col);
		    }
		}
	      */
	      for (int row = 0; row < N_s; row++)
		{
		  for (int col = 0; col < N_N*M_G; col++)
		    {
		      FullSurf(a*N_s + row, col) = SigHat_cont_Surf(row, col);
		    }
		  for (int col = 0; col < N_s; col++)
		    {
		      //negative sign here so I don't need it in DG kernel
		      FullVol(a*N_s + row, col) = (-1.0) * SigHat_cont_Vol(row,col);
		    }
		}
	    } //end of direction 'a' loop

	  //Place FullSurf and FullVol in giant array format.
	  /*
	  for (int row = 0; row < D*N_G; row++)
	    {
	      for (int col = 0; col < N_N*M_G; col++)
		{
		  serial_SigSurf[om*(D*N_G*N_N*M_G) + row*(N_N*M_G) + col] = FullSurf(row,col);
		}
	      for (int col = 0; col < N_s; col++)
		{
		  serial_SigVol[om*(D*N_G*N_s) + row*N_s + col] = FullVol(row,col);
		}
	    }
	  */
	  for (int row = 0; row < D*N_s; row++)
	    {
	      for (int col = 0; col < N_N*M_G; col++)
		{
		  serial_AuxHatSurf[om*(D*N_s*N_N*M_G) + row*(N_N*M_G) + col] = FullSurf(row,col);
		}
	      for (int col = 0; col < N_s; col++)
		{
		  serial_AuxHatVol[om*(D*N_s*N_s) + row*N_s + col] = FullVol(row,col);
		}
	    }
	  delete[] JFace;
	} //end of giant element loop
    }
  //serial_SigSurf and serial_SigVol can now be applied to approximate the gradient
  //over all element interiors. Subroutine concluded
}

void SigFaceMatrices(int M_T, int M_G, int N_s, int N_G, int N_N, int M_s, scalar CHI, fullMatrix<scalar> &phi, fullMatrix<scalar> &dphi, fullMatrix<double> &weight, fullMatrix<scalar> &J, fullMatrix<scalar> &psi, fullMatrix<scalar> &psi_w, fullMatrix<scalar> normals, fullMatrix<scalar> &JF, int* Alt_FaceFromElem, int* BR2_Map, int* invmap, int* Better_InvMap, scalar* invJac, scalar* Minv_global, scalar* PSIxR_Global, scalar* detJ_full, scalar* serial_SigFace_from_DOF)
{
    /*!
    \brief Build operator to populate interface gradient from surrounding DG DOF
    \param[in] M_T total interface count
    \param[in] M_G quadrature nodes per interface
    \param[in] N_s solution points per element
    \param[in] N_G interior quadrature nodes per element
    \param[in] N_N neighbors per element (or faces per element)'
    \param[in] M_s supported basis functions per interface, from one element
    \param[in] CHI controls how much the common interface solution influences SigFace calculation
    \param[in] phi DG basis on reference element
    \param[in] dphi gradient of basis function, ref element, wrt ref coordinates
    \param[in] weight quadrature weights over the reference element
    \param[in] J the determinant of the Jacobian; constant per element
    \param[in] psi trace of DG basis functions on interface
    \param[in] psi_w psi functions multiplied by appropriate face quadrature weights
    \param[in] normals the outward normal from element A along each interface
    \param[in] JF Jacobian of each interface
    \param[in] Alt_FaceFromElem tells each element who its N_N interfaces are
    \param[in] BR2_Map tells each interface who its two elements are (and something about face orientation)
    \param[in] invmap Marc's inverse mapping to relate interfaces to element basis functions
    \param[in] invJac entries of each element's inverse Jacobian. NOT determinant of said inverse.
    \param[in] Minv_global the inverse mass matrix of each element
    \param[in] PSIxR_Global the recovery operator on each interface
    \param[in] detJ_full det(Jacobian) at each quadrature point of each element
    \param[out] serial_SigFace_from_DOF the DOF-to-face gradient operator
  */
   
 
  int verbose = 0; //switch to 1 to get a printout of what is happening in here
  int Special = 0;

  if (verbose == 1) {printf("Entered the SigFaceMatrices subroutine\n");}
  //scalar CHI = 2.0; //the CGR-I penalty parameter
  for (int t = 0; t < M_T; t++)
    {
      int omA = BR2_Map[t*4 + 0];
      int omB = BR2_Map[t*4 + 2];
      int BoundaryTag = 0;
      if (omA == 6 && omB == 32 && Special > 0)
	{
	  verbose = 1;
	  printf("Changed to verbose=1 in GetSigFaceMatrices for omA=6, omB=32 case\n");
	}
      else
	{
	  verbose = 0;
	}

      if (omA == omB)
	{
	  BoundaryTag = 1;
	  if (verbose > 0){printf("interface %d is a boundary interface, will set SigFaceMatrix to zero\n",t);}
	}
      if (verbose == 1){     printf("t=%d, omA=%d, omB=%d\n",t,omA,omB);}
      if (BoundaryTag == 0)
	{
	  //interior interface, take usual approach to build serial_SigFace_from_DOF
	  int HA = -1;//tells the element which of its N_N faces match the global interface t
	  int HB = -1;
	  for (int H = 0; H < N_N; H++)
	    {
	      if (Alt_FaceFromElem[omA*N_N + H] == t)
		{
		  HA = H;
		}
	      if (Alt_FaceFromElem[omB*N_N + H] == t)
		{
		  HB = H;
		}
	    }
	  if (verbose == 1){
	    if (HA == -1)
	      {
		printf("\nCATASTROPHE!!! HA remains -1 in SigFaceMatrices\n");
	      }
	    if (HB == -1)
	      {
		printf("\nHB remains -1 in SigFaceMatrices\n");
	      }
	  }
	  //In periodic case, some interfaces are double-numbered. So, if HA or HB is negative 1,
	  //it means one of the elements does not know that it borders the interface t.

	  if (HB == -1)
	    {
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
				}
			    }
			}
		    }
		}
	    }
	  //      printf("No segfault A\n");
	  //grab the left/right halves of the interface's PSIxR matrix
	  fullMatrix<scalar> PSIxR_omA(M_G, N_s);
	  fullMatrix<scalar> PSIxR_omB(M_G, N_s);
	  for (int g = 0; g < M_G; g++) {
	    for (int col = 0; col < N_s; col++) {
	      PSIxR_omA(g,col) = PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + 0   + col];
	      PSIxR_omB(g,col) = PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + N_s + col]; }}
	  /*
	    printf("PSIxR_omA:\n");
	    for (int row = 0; row < M_G; row++)
	    {
	    printf("Quadrature node %d: ",row);
	    for (int col = 0; col < N_s; col++)
	    {
	    printf("%f, ",PSIxR_omA(row,col));
	    }
	    printf("\n");
	    }
	    printf("PSIxR_omB:\n");
	    for (int row = 0; row < M_G; row++)
	    {
	    printf("Quadrature node %d: ",row);
	    for (int col = 0; col < N_s; col++)
	    {
	    printf("%f, ",PSIxR_omB(row,col));
	    }
	    printf("\n");
	    }
	  */  

	  //grab the inverse mass matrices of the 2 elements
	  fullMatrix<scalar> SSAM_A(N_s, N_s);
	  fullMatrix<scalar> SSAM_B(N_s, N_s);
	  for (int row = 0; row < N_s; row++) {
	    for (int col = 0; col < N_s; col++) {
	      SSAM_A(row,col) = Minv_global[omA*N_s*N_s + col*N_s + row];
	      SSAM_B(row,col) = Minv_global[omB*N_s*N_s + col*N_s + row]; }}
	  /*
	    printf("SSAM_A:\n");
	    for (int row = 0; row < N_s; row++)
	    {
	    printf("row %d: ",row);
	    for (int col = 0; col < N_s; col++)
	    {
	    printf("%f, ",SSAM_A(row,col));
	    }
	    printf("\n");
	    }
	    printf("SSAM_B:\n");
	    for (int row = 0; row < N_s; row++)
	    {
	    printf("row %d: ",row);
	    for (int col = 0; col < N_s; col++)
	    {
	    printf("%f, ",SSAM_B(row,col));
	    }
	    printf("\n");
	    }
	  */

	  //Another issue: Sometimes, invmap doesn't properly pair nodes.
	  //Specifically, for 
	  if (verbose == 1)
	    {
	      printf("The invJac we are using for element B:\n");
	      for (int g = 0; g < N_G; g++)
		{
	      for (int a1 = 0; a1 < D; a1++)
		{
		  for (int a2 = 0; a2 < D; a2++)
		    {
		      printf("g=%d: invJac(a1=%d,a2=%d) = %f\n",g, a1,a2,invJac[omB*N_G*D*D + g*D*D + a1*D + a2]);
		    }
		}
		}
	    }
	  for (int a = 0; a < D; a++)
	    {
	      if (verbose == 1) {printf("interface %d, direction(a)=%d, omA=%d, omB=%d, HA=%d, HB=%d\n", t, a,omA,omB,HA,HB);}
	      //Now, the the gradient contribution
	      //for a particular gradient component

	      //Get the basis derivatives wrt physical coordinate Xa.
	      //invJac contains the necessary information to move between ref and physical coordinates
	      fullMatrix<scalar> dphiA_dXa (N_s, N_G);
	      fullMatrix<scalar> dphiB_dXa (N_s, N_G);
	      for (int k = 0; k < N_s; k++) {
		for (int g = 0; g < N_G; g++) {
		  dphiA_dXa(k,g) = 0.0;
		  dphiB_dXa(k,g) = 0.0;
		  for (int alpha = 0; alpha < D; alpha++) {
		    /*
		    dphiA_dXa(k,g) += dphi(g*D+alpha, k) * invJac[omA*N_G*D*D + g*D*D + a*D + alpha];
		    dphiB_dXa(k,g) += dphi(g*D+alpha, k) * invJac[omB*N_G*D*D + g*D*D + a*D + alpha]; 
		    */
		  dphiA_dXa(k,g) += dphi(g*D+alpha, k) * invJac[omA*N_G*D*D + g*D*D + alpha*D + a];
		  dphiB_dXa(k,g) += dphi(g*D+alpha, k) * invJac[omB*N_G*D*D + g*D*D + alpha*D + a];
		  }}}
		    
	      if (verbose > 0)
		{	
			printf("dphiA_dXa:\n");
			for (int k = 0; k < N_s; k++)
			{
			printf("\tindex %d:\n",k);
			for (int g = 0; g < N_G; g++)
			{
			printf("\t\tnode %d: dphiA_dXa = %f\n", g, dphiA_dXa(k,g));
			}
			}
			printf("dphiB_dXa:\n");
			for (int k = 0; k < N_s; k++)
			{
			printf("\tindex %d\n:",k);
			for (int g = 0; g < N_G; g++)
			{
			printf("\t\tnode %d: dphiB_dXa = %f\n", g, dphiB_dXa(k,g));
			}
			}
	    

			printf("a=%d, no Segfault B\n",a);
		}	  
	      fullMatrix<scalar> subMatrix(M_G, 2*N_s);
	      fullMatrix<scalar> superBlock_A(N_s, 2*N_s);
	      fullMatrix<scalar> superBlock_B(N_s, 2*N_s);
	      fullMatrix<scalar> ssam_by_SuperBlock_A(N_s, 2*N_s);
	      fullMatrix<scalar> ssam_by_SuperBlock_B(N_s, 2*N_s);
	      fullMatrix<scalar> Cont_A(M_G, 2*N_s);
	      fullMatrix<scalar> Cont_B(M_G, 2*N_s);

	      //VolMat_1 x VolMat_2 is a bit like mass matrix, except I take a gradient in direction a
	      fullMatrix<scalar> VolMat_1_A(N_s, N_G);
	      fullMatrix<scalar> VolMat_2_A(N_G, N_s);
	      fullMatrix<scalar> VolMat_1_B(N_s, N_G);
	      fullMatrix<scalar> VolMat_2_B(N_G, N_s);
	      /*
		printf("detJ of omA=%f||detJ of omB=%f\n",J(omA,0), J(omB,0));
		printf("imported phiRef:\n");
		for (int g = 0; g < N_G; g++)
		{
		printf("\tnode %d quad weight = %f:\n",g, weight(g,0));
		for (int k = 0; k < N_s; k++)
		{
		printf("\t\tphi(node=%d,index=%d) = %f\n",g,k,phi(g,k));
		}
		}
	      */
	      for (int k = 0; k < N_s; k++) {
		for (int g = 0; g < N_G; g++) {
		  //VolMat_1_A(k,g) = phi(g,k) * weight(g,0) * J(omA,0);
		  VolMat_1_A(k,g) = phi(g,k) * weight(g,0) * detJ_full[omA*N_G + g];
		  VolMat_2_A(g,k) = dphiA_dXa(k,g);
		  //VolMat_1_B(k,g) = phi(g,k) * weight(g,0) * J(omB,0);
		  VolMat_1_B(k,g) = phi(g,k) * weight(g,0) * detJ_full[omB*N_G + g];
		  VolMat_2_B(g,k) = dphiB_dXa(k,g); }}

	      fullMatrix<scalar> FaceIntegrate_A(N_s, M_G);
	      fullMatrix<scalar> FaceIntegrate_B(N_s, M_G);
	      //start by zeroing these matrices
	      for (int row = 0; row < N_s; row++){
		for (int g = 0; g < M_G; g++){
		  FaceIntegrate_A(row,g) = 0.0;
		  FaceIntegrate_B(row,g) = 0.0; }}
	      if (verbose > 0)
		{
		printf("a=%d, no Segfault C\n",a);
		printf("psi and psi_w structure:\n");
		for (int g = 0; g < M_G; g++)
		{
		for (int m = 0; m < M_s; m++)
		{
		printf("psi(node=%d,index=%d)=%f\t\tpsi_w(node=%d,index=%d) = %f\n", g,m,psi(g,m),g, m, psi_w(g,m));
		}
		}
		printf("face normal in direction a=%d: %f\n", a, normals(a,t));
		printf("Face Jacobian: %f and %f\n", JF(t*2 + 0,0), JF(t*2 + 1,0));
		}
	      //face integration: must fight with supported node issue
	      if (verbose == 1) {printf("FaceIntegrate_A:\n");}
	      for (int m = 0; m < M_s; m++)
		{
		  //int row = invmap[((omA*N_F+0)*M_s*N_N +  HA*M_s + m)*2+0]; //solIDX of the inteface-supported basis function
		  //row = the element-wise node number corresponding to the facial node m, take fc=0
		  int fc = 0;
		  int side = 0; //omA gets side 0
		  int face = ((t*N_F+fc)*2+side)*M_s+m;
		  int Masterf0Index = Better_InvMap[face];
		  int row = Masterf0Index - omA*(N_F*N_s) - fc*N_s; //gives k = element index corresponding to the node
		  if (verbose == 1) {printf("m=%d,row=%d\n",m,row);}
		  for (int g = 0; g < M_G; g++)
		    {
		      FaceIntegrate_A(row,g) = psi_w(g,m) * normals(a,t) * JF(t*2 + 0, 0); } }
	      //Repeat for B element:
	      //	  printf("a=%d, no Segfault D\n",a);
	      if (verbose == 1) {printf("FaceIntegrate_B:\n");}
	      for (int m = 0; m < M_s; m++)
		{
		  //int row = invmap[((omB*N_F+0)*M_s*N_N +  HB*M_s + m)*2+0]; //solIDX of the intefacre-supported basis function
		  int fc = 0;
		  int side = 1; //omB gets side 1
		  int face = ((t*N_F+fc)*2+side)*M_s+m;
		  int Masterf0Index = Better_InvMap[face];
		  int row = Masterf0Index - omB*(N_F*N_s) - fc*N_s; //gives k = element index corresponding to the node
		  if (verbose == 1) {printf("m=%d,row=%d\n",m,row);}
		  for (int g = 0; g < M_G; g++)
		    {
		      //negative on normal because it is outward from A element
		      FaceIntegrate_B(row,g) = psi_w(g,m) * (-normals(a,t)) * JF(t*2 + 1, 0); } }

	      //	  printf("a=%d, no Segfault E\n",a);

	      //grabTrace populates the DG solution along the interface from DG DOF.
	      //Marc's map sturcture:
	      if (verbose > 0)
		{
		  printf("FaceIntegrateA, the whole thing:\n");
		  for (int k = 0; k < N_s; k++)
		    {
		      for (int g = 0; g < M_G; g++)
			{
			  printf("\tFaceIntegrate(k=%d, g=%d) = %f\n", k, g, FaceIntegrate_A(k,g));
			}
		    }
		  printf("FaceIntegrateB, the whole thing:\n");
		  for (int k = 0; k < N_s; k++)
		    {
		      for (int g = 0; g < M_G; g++)
			{
			  printf("\tFaceIntegrate(k=%d, g=%d) = %f\n", k, g, FaceIntegrate_B(k,g));
			}
		    }
		}
	  
	      fullMatrix<scalar> grabTrace_A(M_G, N_s);
	      fullMatrix<scalar> grabTrace_B(M_G, N_s);
	      //zero the grabTrace matrices
	      for (int row = 0; row < M_G; row++)
		{
		  for (int col = 0; col < N_s; col++)
		    {
		      grabTrace_A(row,col) = 0.0;
		      grabTrace_B(row,col) = 0.0;
		    }
		}
	      for (int m = 0; m < M_s; m++) {
		//int row = invmap[((omA*N_F+0)*M_s*N_N +  HA*M_s + m)*2+0]; //solIDX of the intefacre-supported basis function
		int fc = 0;
		int side = 0; //omA gets side 0
		int face = ((t*N_F+fc)*2+side)*M_s+m;
		int Masterf0Index = Better_InvMap[face];
		int row = Masterf0Index - omA*(N_F*N_s) - fc*N_s; //gives k = element index corresponding to the node
		for (int g = 0; g < M_G; g++) {
		  grabTrace_A(g,row) = psi(g,m); }}
	      //repeat for B:
	      for (int m = 0; m < M_s; m++) {
		// int row = invmap[((omB*N_F+0)*M_s*N_N +  HB*M_s + m)*2+0]; //solIDX of the intefacre-supported basis function
		int fc = 0;
		int side = 1; //omB gets side 1
		int face = ((t*N_F+fc)*2+side)*M_s+m;
		int Masterf0Index = Better_InvMap[face];
		int row = Masterf0Index - omB*(N_F*N_s) - fc*N_s; //gives k = element index corresponding to the node
		for (int g = 0; g < M_G; g++) {
		  grabTrace_B(g,row) = psi(g,m); }}
	  
	      //	  printf("a=%d, no Segfault F\n",a);
	  
	      //Assemble the superBlock matrices: 
	      //The easy part is omA dependence on omB and omB depencence on omA.
	      //This is right half of omA superblock and left half of omB superblock
	      fullMatrix<scalar> M_cross_A(N_s, N_s); 
	      fullMatrix<scalar> M_cross_B(N_s, N_s);
	      M_cross_A.gemm(FaceIntegrate_A, PSIxR_omB);
	      M_cross_B.gemm(FaceIntegrate_B, PSIxR_omA);
	      for (int row = 0; row < N_s; row++) {
		for (int col = 0; col < N_s; col++) {
		  superBlock_A(row,col + N_s) = CHI*M_cross_A(row,col); 
		  superBlock_B(row,col + 0  ) = CHI*M_cross_B(row,col);}}
		  
	      //Next: the home element dependence.
	      //The home element contributes a volume integral,
	      //but also sends a trace to the interface
	      fullMatrix<scalar> M_home1_A(N_s, N_s);
	      fullMatrix<scalar> M_home2_A(N_s, N_s);
	      fullMatrix<scalar> M_home1_B(N_s, N_s);
	      fullMatrix<scalar> M_home2_B(N_s, N_s);
	      fullMatrix<scalar> M_faceDiff_A(M_G, N_s);
	      fullMatrix<scalar> M_faceDiff_B(M_G, N_s);
	      for (int row = 0; row < M_G; row++) {
		for (int col = 0; col < N_s; col++) {
		  //faceDiff accounts for recovered solution minus home trace.
		  //That difference must be integrated subject to FaceIntegrate
		  M_faceDiff_A(row,col) = PSIxR_omA(row,col) - grabTrace_A(row,col);
		  M_faceDiff_B(row,col) = PSIxR_omB(row,col) - grabTrace_B(row,col); }}
	  
	      M_home1_A.gemm(VolMat_1_A, VolMat_2_A);
	      M_home2_A.gemm(FaceIntegrate_A, M_faceDiff_A);
	      M_home1_B.gemm(VolMat_1_B, VolMat_2_B);
	      M_home2_B.gemm(FaceIntegrate_B, M_faceDiff_B);

	      //Now, relay the home-element dependence to superBlock
	      for (int row = 0; row < N_s; row++) {
		for (int col = 0; col < N_s; col++) {
		  superBlock_A(row, col + 0  ) = M_home1_A(row,col) + CHI*M_home2_A(row,col);
		  superBlock_B(row, col + N_s) = M_home1_B(row,col) + CHI*M_home2_B(row,col); }}
	  

	      //Multiply by inverse mass matrix in each element
	      ssam_by_SuperBlock_A.gemm(SSAM_A, superBlock_A);
	      ssam_by_SuperBlock_B.gemm(SSAM_B, superBlock_B);
	  
	      //Multiply by PSIxR operator in each element
	      Cont_A.gemm(PSIxR_omA, ssam_by_SuperBlock_A);
	      Cont_B.gemm(PSIxR_omB, ssam_by_SuperBlock_B);
	  
	      //Add the contributions together to get the (t,a) submatrix
	      for (int row = 0; row < M_G; row++) {
		for (int col = 0; col < 2*N_s; col++) {
		  subMatrix(row,col) = Cont_A(row,col) + Cont_B(row,col); }}


	      if (verbose > 0)
		{
		  printf("Testing the SigFace operation for t=%d, a=%d\n",t,a);
		  scalar UA[N_s];
		  scalar UB[N_s];
		  for (int k = 0; k < N_s; k++)
		    {
		      UA[k] = 1.0;
		      UB[k] = 0.0;
		    }
		  for (int g = 0; g < M_G; g++)
		    {
		      scalar sum = 0;
		      for (int k = 0; k < N_s; k++)
			{
			  sum += subMatrix(g, k) * UA[k];
			  sum += subMatrix(g, k) * UB[k];
			}
		      printf("g=%d: derivative=%f\n", g,sum);
		    }
		}
	  
	      //	  printf("a=%d, no Segfault G\n",a);
	      //relay the output to global serial form:
	      for (int g = 0; g < M_G; g++) {
		for (int col = 0; col < 2*N_s; col++) {
		  serial_SigFace_from_DOF[t*D*M_G*2*N_s + a*M_G*2*N_s + g*2*N_s + col] = subMatrix(g,col); }}
	  
	      // printf("a=%d, no Segfault H\n",a);
	    } //end a<D component loop
	}
      else
	{
	  if (verbose > 0){printf("Boundary interface, settig SigFace matrix to zero\n");}
	  for (int a = 0; a < D; a++){
	    for (int g = 0; g < M_G; g++){
	      for (int col = 0; col < 2*N_s; col++){
		serial_SigFace_from_DOF[t*D*M_G*2*N_s + a*M_G*2*N_s + g*2*N_s + col] = 0.0;}}}
	}
		
    }  //end t < M_T loop
} //end subroutine

/*
void SigSurfOneElementOneFace(int tGlo, int om, int M_s, int* invmap, fullMatrix<scalar> &OUT)
{

  //Supplement to calculate integral over (el e, local face l) [All phi * tilde{U}dot Normal_component ds]
  //The information sent to this routine is specific to one interface of one element.
  //Output must be assembled in a parent routine
  fullMatrix<scalar> ROG (N_s, M_s);
  int solIndex[M_s];
  for (int m = 0; m < M_s; m++)
    {
      //where I have zero, that's where the field index fc usually goes
      solIndex[m] = invmap[((om*N_F+0)*M_s*N_N+i)*2+0];
    }
  //Output matrix is ROG by PSIwPlus
  OUT.GEMM(ROG, PSIwPlus);
  //The output of this routine is the matrix OUT, subroutine concluded
}
*/
