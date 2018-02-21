/*!
  \file dg_functions.cc
  \brief Functions definitions for the DG method
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include "dg_functions.h"
void dg_jac_elements_fast(const int N_G, const int N_E, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &dphi, fullMatrix<scalar> &J){
  /*!
    \brief Get the Jacobian of the elements (only, which is why it's "fast"). Also fast because it stores only one detJ value for the entire element (using full Jacobian matrix at first quadrature point), which is a questionable practice, but we can usually get away with it.
    \param[in] N_G number of gaussian nodes per element
    \param[in] N_E number of elements
    \param[in] XYZNodes coordinates of element nodes
    \param[in] dphi derivative of the polynomial basis
    \param[out] J Jacobian of each element
  */

  fullMatrix<scalar> Jac(N_E*D,N_G*D);
  Jac.gemm(XYZNodes.transpose(),dphi.transpose());
  scalar det = 0.0;
  
  for(int e = 0; e < N_E; e++){
#ifdef ONED
    det = sqrt(Jac(e,0)*Jac(e,0));
#elif TWOD
    det = Jac(e*D+0,0)*Jac(e*D+1,1)-Jac(e*D+1,0)*Jac(e*D+0,1);
#endif
#ifdef THREED
    //Use Cramer's rule to get the determinant (at least, I think it's called Cramer's)
    /*
      det[a b c; 
          d e f; 
	  g h i] = a(ei-fh) - b(di-fg) + c(dh-eg)
     */
    scalar m11 = Jac(e*D+0,0); scalar m12 = Jac(e*D+0,1); scalar m13 = Jac(e*D+0,2);
    scalar m21 = Jac(e*D+1,0); scalar m22 = Jac(e*D+1,1); scalar m23 = Jac(e*D+1,2);
    scalar m31 = Jac(e*D+2,0); scalar m32 = Jac(e*D+2,1); scalar m33 = Jac(e*D+2,2);
    det = m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 - m22*m31);
    /*
    int c1 = 0;
    int c2 = 1;
    int c3 = 2;
    int r1 = 0;
    int r2 = 1;
    int r3 = 2;
    det  = Jac(e*D + r1,c1) * (Jac(e*D + r2,c2)*Jac(e*D + r3,c3)-Jac(e*D + r2,c3)*Jac(e*D + r2,c2));
    det -= Jac(e*D + r1,c2) * (Jac(e*D + r2,c1)*Jac(e*D + r3,c3)-Jac(e*D + r2,c3)*Jac(e*D + r3,c1));
    det += Jac(e*D + r1,c3) * (Jac(e*D + r2,c1)*Jac(e*D + r3,c2)-Jac(e*D + r2,c2)*Jac(e*D + r3,c1));
    */
#endif
    J(e,0) =  fabs(det);
  }
}

void dg_detJ_OneElement(const int N_G, scalar* xGeo, fullMatrix<scalar> &dphi, scalar* detJ)
{
  /*!
    \brief get the determinant of jacobian at quadrature points for a given element
    \param[in] N_G quadrature points per element (doesn't have to match with main.cc)
    \param[in] xGeo physical coordinates of the element's nodes; Ns*D 
    \param[in] dphi the gradient matrix of solution basis wrt reference coordinates; (NG*D) x Ns
    \param[out] detJ the determinant of Jacobain; this routine calculates N_G values
  */
  int verbose = 0;
  int N_s = dphi.size2();
  if (verbose > 0){printf("In dg_detJ_OneElement: N_s=%d\n",N_s);}
  fullMatrix<scalar> Jac; //DxD matrix at a given quadrature node
  Jac.resize(D,D);
  for (int g = 0; g < N_G; g++)
    {
      for (int row = 0; row < D; row++)
	{
	  for (int col = 0; col < D; col++)
	    {
	      Jac(row,col) = 0.0;
	      for (int k = 0; k < N_s; k++)
		{
		  Jac(row,col) += xGeo[k*D + row] * dphi(g*D+col, k);
		  //	  printf("xGeo(k=%d, a=%d)=%f, dphi(g=%d, k=%d, a=%d)=%f
		}
	    }
	}
      //We now have the geometric Jacobian AT 1 QUADRATURE POINT.
      //Now, take determinant
      scalar detJ_local;
      switch(D)
	{
	case 1: {detJ_local = fabs(Jac(0,0)); break;}
	case 2:
	  {
	    detJ_local = fabs(Jac(0,0)*Jac(1,1) - Jac(0,1)*Jac(1,0));
	    break;
	  }
	case 3:
	  {
	    detJ_local =  Jac(0,0) * (Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1));
	    detJ_local -= Jac(0,1) * (Jac(1,0)*Jac(2,2) - Jac(1,2)*Jac(2,0));
	    detJ_local += Jac(0,2) * (Jac(1,0)*Jac(2,1) - Jac(1,1)*Jac(2,0));
	    detJ_local = fabs(detJ_local);
	    break;
	  }
	default:
	  {
	    detJ_local = 0.0;
	    break;
	  }
	}
      detJ[g] = detJ_local;
      if (verbose > 0){printf("In dg_detJ_OneElement: detJ[g=%d]=%f\n", g, detJ[g]);}
    }
}


void dg_jacobians_elements(const int N_G, const int N_E, const int Ne_AUG, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &XYZNodes_extended, fullMatrix<scalar> &dphi, fullMatrix<scalar> &Jac, fullMatrix<scalar> &invJac, fullMatrix<scalar> &J, fullMatrix<scalar> &invJ){
  /*!
    \brief Get the Jacobian and inverse Jacobian of the elements, ghost and flesh
    \param[in] N_G number of gaussian nodes per element
    \param[in] N_E number of elements
    \param[in] Ne_AUG ghost+flesh elements
    \param[in] XYZNodes coordinates of element nodes
    \param[in] XYZNodes_extended the coordinates of element nodes, flesh+ghost elemtns
    \param[in] dphi derivative of the polynomial basis
    \param[out] Jac Jacobian of each element
    \param[out] invJac inverse Jacobian of each element
    \param[out] J Jacobian of each element
    \param[out] invJ inverse Jacobian of each element
  */
  int verbose = 0;
  int Special = 0;
  //Jac.gemm(XYZNodes.transpose(),dphi.transpose());
  Jac.gemm(XYZNodes_extended.transpose(),dphi.transpose());
  scalar det = 0.0;
  
  //for(int e = 0; e < N_E; e++)
  for(int e = 0; e < Ne_AUG; e++)
    {
      if (e == 32 && Special > 0)
	{
	  verbose = 1;
	  printf("In dg_jacobians_elements: set verbose=1 for om==32\n");
	}
      else
	{
	  verbose = 0;
	}
      if (verbose == 1)
	{
	  printf("Local element's XYZNodes matrix\n");
	  for (int k = 0; k < XYZNodes_extended.size1(); k++)
	    {
	      printf("k=%d: x,y,x = (%f, %f, %f)\n", k, XYZNodes_extended(k, e*D+0),XYZNodes_extended(k, e*D+1),XYZNodes_extended(k, e*D+2));
	    }
	}

#ifdef ONED
    det = sqrt(Jac(e,0)*Jac(e,0));   
#elif TWOD
    det = Jac(e*D+0,0)*Jac(e*D+1,1)-Jac(e*D+1,0)*Jac(e*D+0,1);
#elif THREED
    //Use Cramer's rule to get the determinant (at least, I think it's called Cramer's)
    /*
      det[a b c; 
          d e f; 
	  g h i] = a(ei-fh) - b(di-fg) + c(dh-eg)
     */
    scalar m11 = Jac(e*D+0,0); scalar m12 = Jac(e*D+0,1); scalar m13 = Jac(e*D+0,2);
    scalar m21 = Jac(e*D+1,0); scalar m22 = Jac(e*D+1,1); scalar m23 = Jac(e*D+1,2);
    scalar m31 = Jac(e*D+2,0); scalar m32 = Jac(e*D+2,1); scalar m33 = Jac(e*D+2,2);
    det = m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 - m22*m31);
#endif
    J(e,0)    =  fabs(det); //shitty, element-constant Jacobian
    invJ(e,0) =  fabs(1.0/det);
    for(int g = 0; g < N_G; g++){
      //For effective invJac calculation at quadrature nodes,
      //must use the actual Jacoobian at that node 
      //and not assume constant detJ over entire element
#ifdef ONED
      scalar detLocal = fabs(Jac(e,g));
      invJac(e,g) = 1.0/Jac(e,g);
#elif TWOD
      //no messing around here; perform true inversion
      fullMatrix<double> JacNode(D,D);
      fullMatrix<double> invJacNode(D,D);
      for (int a1 = 0; a1 < D; a1++)
	{
	  for (int a2 = 0; a2 < D; a2++)
	    {
	      //Not sure if this indexing is right
	      JacNode(a1,a2) = Jac(e*D + a1, g*D + a2);
	    }
	}
      //Invert JacNode, store result in invJacNode
      //printf("e=%d(N_E=%d, Ne_AUG=%d), g=%d, preparing to invert JacNode\n",e,N_E, Ne_AUG,g);
      
      JacNode.invert(invJacNode);
      /*
      if (fabs(JacNode(1,1)) < pow(10,-12))
	{
	  printf("Shitty matrix, terminating invJac loop early in dg_jacobians_elements\n");
	  e=Ne_AUG;
	  g = N_G;
	  printf("Here's the shitty matrix:\n");
	  for (int a1 = 0; a1 < D; a1++){
	    printf("row %d:  ",a1);
	    for (int a2 = 0; a2 < D; a2++){
	      printf("%f, ",JacNode(a1,a2));}
	    printf("\n");
	  }
	}
      */
      //Now, relay the result to global storage,
      //using same jumbled ordering as 2D case
      for (int a1 = 0; a1 < D; a1++)
	{
	  for (int a2 = 0; a2 < D; a2++)
	    {
	      //not sure if this indexing is right
	      //invJac(e*D + a2, g*D + a1) = invJacNode(a1,a2);
	      invJac(e*D + a2, g*D + a1) = invJacNode(a2,a1); //Correct, verified on non-uniform quad mesh
	    }
	}
      /*
	/Original approach: use element-constant det for inversion process
      invJac(e*D+0,g*D+0) =  1.0/det*Jac(e*D+1,g*D+1);
      //invJac(e*D+1,g*D+0) = -1.0/det*Jac(e*D+0,g*D+1); //caused problem
      //invJac(e*D+0,g*D+1) = -1.0/det*Jac(e*D+1,g*D+0); //caused problem
      invJac(e*D+1,g*D+0) = -1.0/det*Jac(e*D+1,g*D+0); //original
      invJac(e*D+0,g*D+1) = -1.0/det*Jac(e*D+0,g*D+1); //original
      invJac(e*D+1,g*D+1) =  1.0/det*Jac(e*D+0,g*D+0);
      */
#elif THREED
      fullMatrix<double> JacNode(D,D);
      fullMatrix<double> invJacNode(D,D);
      for (int a1 = 0; a1 < D; a1++)
	{
	  for (int a2 = 0; a2 < D; a2++)
	    {
	      //Not sure if this indexing is right
	      JacNode(a1,a2) = Jac(e*D + a1, g*D + a2);
	 
	      if (verbose>0){printf("at g=%d: Jac((%d,%d) = %f\n", g, a1,a2,JacNode(a1,a2));}
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
	      //invJac(e*D + a2, g*D + a1) = invJacNode(a1,a2);
	      invJac(e*D + a2, g*D + a1) = invJacNode(a2,a1);
	      if (verbose>0){printf("at g=%d: invJac((%d,%d) = %f\n", g, a2,a1,invJacNode(a2,a1));}
	    }
	}
#endif
    }
  }

  //Jacobian and inverse of the jacobian have now been populated
  //at all quadrature points in all elements
  if (verbose == 1)
    {
      //The Jacobian of the first element:
      printf("---Jacobian of first element:---\n");
      for (int g = 0; g < N_G; g++)
	{
	  printf("\tQuadrature node %d:\n",g);
	  for (int row = 0; row < D; row++)
	    {
	      for (int col = 0; col < D; col++)
		{
		  printf("\t\td(x%d)_d(xi%d)=%f,\t",row,col,Jac(0*D+row, g*D+col));
		}
	      printf("\n");
	    }
	}
      
    }
}


void dg_jacobians_face(const int M_T, fullMatrix<scalar> &XYZNodesF, fullMatrix<scalar> &dpsi, const fullMatrix<scalar> &normals, fullMatrix<scalar> &JacF, fullMatrix<scalar> &JF, fullMatrix<scalar> &invJF){

  /*!
    \brief Get the Jacobian and inverse Jacobian of the interfaces
    \param[in] M_T number of elements
    \param[in] XYZNodesF coordinates of element nodes
    \param[in] dpsi derivative of the interface polynomial basis
    \param[in] normals matrix of face normals.
    \param[out] JacF Jacobian of each interface
    \param[out] JF Jacobian of each interface
    \param[out] invJF inverse Jacobian of each interface
  */
  //PEJ editing 09/15/2017 to account for 3D case
  int verbose = 0;
  JacF.gemm(XYZNodesF.transpose(),dpsi.transpose());
  for(int t = 0; t < M_T; t++){
    for(int d = 0; d < 2; d++){
      scalar det = 0.0;
#ifdef ONED
      det = 1; // in 1D set JF to 1
       JF(t*2+d,0)    =  sqrt(det);
       invJF(t*2+d,0) =  1.0/sqrt(det);
#elif TWOD
      for(int alpha = 0; alpha < D; alpha ++){ det += JacF((t*2+d)*D+alpha,0)* JacF((t*2+d)*D+alpha,0);}
       JF(t*2+d,0)    =  sqrt(det);
       invJF(t*2+d,0) =  1.0/sqrt(det);
#elif THREED
       //PEJ 09/22/2017: Build face jacobian based on node locations, similar
       //to building element jacobians.
       //dpsi(g*DF + alpha, node)
       //XYZNodes(node,(t*2+side)*D + alpha)
       scalar sum = 0;
       scalar LocalJac[3][3]; //3x3 jacobian whose determinant will be used to calculate detJacF.
       //xi,eta = face reference coordinates, part of psi definition.
       //Also, we are only accessing the g=0 portion (first two columns) of JacF under the assumption
       //that the jacobian is uniform across the face
       LocalJac[0][0] = JacF((t*2 + d)*D + 0 , 0); //dx_dxi
       LocalJac[0][1] = JacF((t*2 + d)*D + 0 , 1); //dx_deta
       LocalJac[0][2] = normals(0, t); //face normal x direction
       
       LocalJac[1][0] = JacF((t*2 + d)*D + 1 , 0); //dy_dxi
       LocalJac[1][1] = JacF((t*2 + d)*D + 1 , 1); //dy_deta
       LocalJac[1][2] = normals(1, t); //face normal y direction
       
       LocalJac[2][0] = JacF((t*2 + d)*D + 2 , 0); //dz_dxi
       LocalJac[2][1] = JacF((t*2 + d)*D + 2 , 1); //dz_deta
       LocalJac[2][2] = normals(2, t); //face normal z direction
       //Now, we take the determinant of this matrix 
       //to assist with integration over the face wherever necessary.
       //I'm not sure what the inverse of this Jacobian indicates.
       
       JF(t*2+d,0)  = LocalJac[0][0] * (LocalJac[1][1]*LocalJac[2][2] - LocalJac[1][2]*LocalJac[2][1]);
       JF(t*2+d,0) -= LocalJac[0][1] * (LocalJac[1][0]*LocalJac[2][2] - LocalJac[1][2]*LocalJac[2][0]);
       JF(t*2+d,0) += LocalJac[0][2] * (LocalJac[1][0]*LocalJac[2][1] - LocalJac[1][1]*LocalJac[2][0]);
       JF(t*2+d,0) = fabs(JF(t*2+d,0));
       //the determinant of a matrix is the same as the inverse of the determinant of its inverse
       invJF(t*2+d,0) = 1.0 / JF(t*2+d,0);
       /*       
       //This part is rather confusing; physical coordinates are 3D but face's reference
       //coordinates are 2D, so 3x3 Jacobian cannot be formed at each point. Will blindly use
       //Marc's formula above
       for(int alpha = 0; alpha < D; alpha ++){ det += JacF((t*2+d)*D+alpha,0)* JacF((t*2+d)*D+alpha,0);}
       //originally took square root here, but just taking fabs gives good results for cube elements
       JF(t*2+d,0) = fabs(det);
       invJF(t*2+d,0) = 1.0 / fabs(det);
       */
#endif

    }
  }
  //Each face now has a jacobian (good for integrating surface terms( and inv Jacobian (don't know what that is for)
  if (verbose == 1)
    {
      for (int t = 0; t < M_T; t++)
	{
	  printf("Interface detJ, interface number %d: JF=%f, invJF=%f\n",t, JF(t*2+0,0),invJF(t*2+0,0));
	}
    }
}


void dg_inverse_mass_matrix(const int order, const int elem_type, const std::string getElemType, const int N_s, const int N_E, const int Ne_AUG, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &XYZNodes_extended, int QuadRule, scalar* detJ_full, scalar* Minv){
  /*!
    \brief Get the inverse mass matrix for each element
    \param[in] order DG order
    \param[in] elem_type type of element (numerical value)
    \param[in] getElemType type of element (string)
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] Ne_AUG, flesh+ghost elements
    \param[in] XYZNodes coordinates of element nodes
    \param[in] XYZNodes_extended coordinates of element nodes, flesh and ghost elements
    \param[in] detJ_full the determinant of geometric Jacobian at each quadrature point, each element
    \param[out] Minv inverse mass matrix for each element
  */
  int verbose = 0;
  //PEJ: argument sent to gaussIntegration should be more than order*2+1
  //if we are dealing with deformed elements. In fact, I think it would
  //be a good idea to use ridiculously high order quadrature here in future,
  //just to be safe.

  const polynomialBasis *basis  = polynomialBases::find (elem_type);  // for the element
  fullMatrix<double> points, weight;
  if     (getElemType == "tri") gaussIntegration::getTriangle(QuadRule, points, weight); //order*2+1 because you are integrating phi*phi
  else if(getElemType == "qua") gaussIntegration::getQuad(QuadRule, points, weight);
  else if(getElemType == "lin") gaussIntegration::getLine(QuadRule, points, weight);
  else if(getElemType == "hex") gaussIntegration::getHexahedron(QuadRule, points, weight);

  // Number of integration points           (g index)
  int N_G = points.size1();    

  // Get basis functions
  fullMatrix<scalar> phi (N_G,N_s);
  fullMatrix<double> phiD (N_G,N_s);
  fullMatrix<scalar> dphi(N_G*D,N_s);
  basis->f (points, phiD);
  for(int g = 0; g<N_G; g++){
    for(int i = 0; i < N_s; i++){
      phi(g,i) = (scalar)phiD(g,i);
    }	  
  }
  fullMatrix<scalar> phi_w (N_G,N_s);
  for(int g = 0; g < N_G; g++){
    for(int i = 0; i < N_s; i++){
      phi_w(g,i) = (scalar)phi(g,i) * weight(g,0);
    }
  }

  double grads[N_s][3];  
  for(int g = 0; g<N_G; g++){
    basis->df(points(g,0),points(g,1),points(g,2),grads);
    for(int alpha = 0; alpha < D; alpha ++){
      for(int i = 0; i < N_s; i++){
	dphi(g*D+alpha,i) = (scalar)grads[i][alpha];  // see paper for indexing p.6
      }	  
    }    
  }

  //PEJ 10/31/2017: I don't know if it matters,
  //but to put my mind at ease, I'm introducing a 
  //normalization (norm) to keep matrix values around
  //unity before inversion, then I'll make a post-inversion correction
  
  // Calculate the jacobians
  //To achieve more accurate mass matrix calculation in future,
  //we could instead store the Jacobian at each quadrature point,
  //using dg_jacobians, not the fast option.
  //fullMatrix<scalar> J(N_E,1);            // determinant of the Jacobian
  fullMatrix<scalar> J(Ne_AUG,1);            // determinant of the Jacobian
  //dg_jac_elements_fast(N_G, N_E, XYZNodes, dphi, J);
  dg_jac_elements_fast(N_G, Ne_AUG, XYZNodes_extended, dphi, J);

  //for(int e = 0; e < N_E; e++)
  for(int e = 0; e < Ne_AUG; e++)
    {
      scalar norm = 0.0;
      for (int g = 0; g < N_G; g++)
	{
	  norm += fabs(detJ_full[e*N_G+g]);
	}
      //for small elements detJ is very small; it's average detJ value over element
      norm = norm / (N_G+0.0);
      // Init de psi_w_Jac qui va etre construit pour chaque element
      fullMatrix<scalar> phi_w_Jac = phi_w;
      // Construction de la matrice de masse: M(i,j) = \int psi_i psi_j
      fullMatrix<scalar> M(N_s,N_s);
      scalar jac = J(e,0);
      for(int g = 0; g < N_G; g++){
	for(int i = 0; i < N_s; i++){
	  //phi_w_Jac(g,i) *= jac;
	  phi_w_Jac(g,i) *= detJ_full[e*N_G + g] / norm; //account for element-wise variation
	}
      }
      M.gemm(phi.transpose(),phi_w_Jac); //this is the element's mass matrix
      if (verbose == 1)
	{
	  printf("Element %d Mass Matrix, divided by norm=%f:\n",e,norm);
	  for (int row = 0; row < N_s; row++)
	    {
	      printf("row %d: ",row);
	      for (int col = 0; col < N_s; col++)
		{
		  printf("%f, ",M(row,col));
		}
	      printf("\n");
	    }
	}
    
      // Inverser la matrice de masse.
      fullMatrix<scalar> M_inv(N_s,N_s);
      M.invert(M_inv);
      if (verbose == 1)
	{
	  printf("Element %d INVERSE of Mass Matrix, multiplied by norm=%f:\n",e,norm);
	  for (int row = 0; row < N_s; row++)
	    {
	      printf("row %d: ",row);
	      for (int col = 0; col < N_s; col++)
		{
		  printf("%f, ",M_inv(row,col));
		}
	      printf("\n");
	    }
	  printf("Element %d INVERSE of TRUE Mass Matrix:\n",e);
	  for (int row = 0; row < N_s; row++)
	    {
	      printf("row %d: ",row);
	      for (int col = 0; col < N_s; col++)
		{
		  printf("%f, ",M_inv(row,col) / norm);
		}
	      printf("\n");
	    }
	}
      for(int k = 0; k < N_s; k++){
	for(int i = 0; i < N_s; i++){
	  //	Minv[(e*N_s+k)*N_s+i] = M_inv(i,k);
	  Minv[(e*N_s+k)*N_s+i] = M_inv(i,k) / norm;
	}
      } //k loop
    } //element loop to Ne_AUG
}


int LocalSeekGhost(const std::vector<int> &ghostElement_Key, const std::vector<int> &ghostElement_Val, const std::vector<int> &ghostElement_HNe, int GmshId_on, int GmshId_off)
{
  int ghostId = -1;
  for (int j = 0; j < ghostElement_Key.size(); j++)
    {
      int key = ghostElement_Key[j];
      if (key == GmshId_off)
	{
	  //We have found a key value corresponding to the flesh element with gmsh address GmshId_off.
	  //Problem: there may be multiple entries in _ghostElement_Key that match GmshId_off (partition corner intersections).
	  //Secondary check to see that elId=GmshId_on and elIndex=_ghostElement_HNe[j] are same element.
	  int om_off = ghostElement_Val[j]; //myid index of the ghost replica of GmshId_off, corresponds to _ghostElement_Key[j]
	  int GmshId_HNe = ghostElement_HNe[j];
	  if (GmshId_on == GmshId_HNe)
	    {
	      ghostId = om_off;
	    }
	  else
	    {
	      //NO MATCH! try again :)
	    }
	 
	}
    }
  return ghostId;
}

void dg_mappings(const int myid, const int M_s, const int M_T, const int N_s, const int N_E, const int N_N, const std::vector<simpleInterface> &interfaces, const std::map<int,int> &ElementMap, const std::vector<int> &ghostElement_Key, const std::vector<int> &ghostElement_Val, const std::vector<int> &ghostElement_HNe, const std::vector<std::vector<int> > &closures, int* map, int* invmap, int* Alt_FaceFromElem_serial){
  /*!
    \brief Map interfaces and elements to each other
    \param[in] myid my processor ID (eg. MPI)
    \param[in] M_s number of nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] N_N number of neighbors per element
    \param[in] interfaces vector of all the interfaces
    \param[in] ElementMap map of element ID to element index
    \param[in] ghostElement_Key see simpleMesh.cc
    \param[in] ghostElement_Val see simpleMesh.cc
    \param[in] ghostElement_HNe
    \param[in] closures closures
    \param[out] map map from elements to interfaces
    \param[out] invmap map from interfaces to elements
    \section Description
    Map the element matrix to the interfaces matrix and
    vice-versa. This is a modified version of what I had
    previously. I think it works better and is more efficient.
  */

  //PEJ 05/23/2017: Modifying this function to also output FaceFromElem,
  //which will give the global interface address
  //corresponding to one of each element's N_N interfaces

  //PEJ 11/17/2017: Mapping is flawed for HiOCFD5 vortex meshes.
  //I will fix this issue now.
  

  // Initialize a counter for the inverse map
  int* cnt = new int[N_F*N_E]; for(int k=0;k<N_E*N_F;k++){cnt[k]=0;}

  // Loop over all the interfaces
  for(int t = 0; t < M_T; t++){
    // printf("Inside dg_mappings: face=%d\n",t);

    const simpleInterface &face = interfaces[t];   // get the interface
    const simpleElement *el1 = face.getElement(0); // get the element to the left
    const simpleElement *el2 = face.getElement(1); // get the element to the right
    int e1 = ElementMap.at(el1->getId());
    int el2num = 0; // do this like in buildNeighbors
    if(el2->getPartition()==myid){ el2num = ElementMap.at(el2->getId());}
    else                         {el2num = LocalSeekGhost(ghostElement_Key, ghostElement_Val, ghostElement_HNe, el1->getId(), el2->getId());} // refer to a back column of U
    //else                         { el2num = ghostElementMap.at(el2->getId());} // refer to a back column of U
    int id1 = face.getClosureId(0);
    int id2 = face.getClosureId(1);
    const std::vector<int> &cl1 = closures[id1];
    const std::vector<int> &cl2 = closures[id2];
    /*
    printf("id1=%d, id2=%d\n",id1,id2);
    printf("closure::\n");
    for (int m = 0; m < M_s; m++)
      {
	printf("cl1[%d]=%d     cl2[%d]=%d\n",m,cl1[m],m,cl2[m]);
      }
    */
    for(int fc = 0; fc < N_F; fc++){
      for(int j = 0; j < M_s; j++){

	// Map from U to UF
	map[((t*N_F+fc)*2+0)*M_s+j] = (e1*N_F+fc)*N_s+cl1[j];
	map[((t*N_F+fc)*2+1)*M_s+j] = (el2num*N_F+fc)*N_s+cl2[j];

	// Inverse map
	invmap[((e1*N_F+fc)*M_s*N_N+cnt[e1*N_F+fc])*2+0] = cl1[j];
	invmap[((e1*N_F+fc)*M_s*N_N+cnt[e1*N_F+fc])*2+1] = ((t*N_F+fc)*2+0)*M_s+j;

	//printf("f=%d,j=%d: just set invmap[[(%d*N_F+%d)*M_s*N_N + %d)*2+0] = %d\n", fc,j,e1,fc,cnt[e1*N_F+fc], cl1[j]);
	//printf("f=%d,j=%d: just set invmap[[(%d*N_F+%d)*M_s*N_N + %d)*2+1] = %d\n", fc,j,e1,fc,cnt[e1*N_F+fc], ((t*N_F+fc)*2+0)*M_s+j);
	//PEJ: I think id is always less than N_N for first element (A), geq N_N for second element (B)
	/*
	Alt_FaceFromElem_serial[e1*N_N + id1%N_N] = t;
	printf("Set Alt_FaceFromElem_serial[%d*N_N + %d] = %d\n",e1,id1%N_N, t);
	*/
	if (fc == 0)
	  {
	    //	    Alt_FaceFromElem_serial[e1*N_N + id1%N_N] = t;
	    //printf("Set Alt_FaceFromElem_serial[%d*N_N + %d] = %d\n",e1,id1%N_N, t);
	    Alt_FaceFromElem_serial[e1*N_N + cnt[e1*N_F+fc]/M_s] = t;
	    //	   printf("Set Alt_FaceFromElem_serial[%d*N_N + %d] = %d\n",e1,cnt[e1*N_F+fc]/M_s, t);
	  }
	cnt[e1*N_F+fc]++;
	// If this face is not on the boundary (belongs to my
	// partition and physical==0), then you can use the interface
	// to map to an element on my partition	
	if((face.getPhysicalTag()==0)&&(el2->getPartition()==myid)){
	  int e2 = ElementMap.at(el2->getId());
	  invmap[((e2*N_F+fc)*M_s*N_N+cnt[e2*N_F+fc])*2+0] = cl2[j];
	  invmap[((e2*N_F+fc)*M_s*N_N+cnt[e2*N_F+fc])*2+1] = ((t*N_F+fc)*2+1)*M_s+j;
	  //printf("f=%d,j=%d: just set invmap[[(%d*N_F+%d)*M_s*N_N + %d)*2+0] = %d\n", fc,j,e2,fc,cnt[e2*N_F+fc], cl2[j]);
	  //printf("f=%d,j=%d: just set invmap[[(%d*N_F+%d)*M_s*N_N + %d)*2+1] = %d\n", fc,j,e2,fc,cnt[e2*N_F+fc], ((t*N_F+fc)*2+1)*M_s+j);
	  if (fc == 0)
	    {
	      Alt_FaceFromElem_serial[e2*N_N + cnt[e2*N_F+fc]/M_s] = t;
	  //  printf("Set Alt_FaceFromElem_serial[%d*N_N + %d] = %d\n",e2,cnt[e2*N_F+fc]/M_s, t);
	    }
	  cnt[e2*N_F+fc]++;
	  /*
	    Alt_FaceFromElem_serial[e2*N_N + id2%N_N] = t;
	  printf("Set Alt_FaceFromElem_serial[%d*N_N + %d] = %d\n",e2,id2%N_N, t);
	  */
	}
      }
    }
  }
  delete[] cnt;
}

void dg_mappings_PERIFIX(const int myid, const int M_s, const int M_T, const int N_s, const int N_E, const int N_N, const std::vector<simpleInterface> &interfaces, const std::map<int,int> &ElementMap, const std::map<int,int> &ghostElementMap, const std::vector<std::vector<int> > &closures, int* map, int* invmap, int* Alt_FaceFromElem_serial){
  /*!
    \brief Map interfaces and elements to each other
    \param[in] myid my processor ID (eg. MPI)
    \param[in] M_s number of nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] N_N number of neighbors per element
    \param[in] interfaces vector of all the interfaces
    \param[in] ElementMap map of element ID to element index
    \param[in] ghostElementMap map of element ID to element index for elements in other partitions
    \param[in] closures closures
    \param[out] map map from elements to interfaces
    \param[out] invmap map from interfaces to elements
    \section Description
    Map the element matrix to the interfaces matrix and
    vice-versa. This is a modified version of what I had
    previously. I think it works better and is more efficient.
  */

  //PEJ 05/23/2017: Modifying this function to also output FaceFromElem,
  //which will give the global interface address
  //corresponding to one of each element's N_N interfaces

  //New goal: build map/invmap in such a way that the duplicate
  //interfaces apprearing for periodic boundary elements
  //appear at end of each element's invmap collection

  int verbose = 1;

  // Initialize a counter for the inverse map
  int* cnt = new int[N_F*N_E]; for(int k=0;k<N_E*N_F;k++){cnt[k]=0;}

  // Loop over all the interfaces
  for(int t = 0; t < M_T; t++){
    // printf("Inside dg_mappings: face=%d\n",t);

    const simpleInterface &face = interfaces[t];   // get the interface
    const simpleElement *el1 = face.getElement(0); // get the element to the left
    const simpleElement *el2 = face.getElement(1); // get the element to the right
    int e1 = ElementMap.at(el1->getId());
    int el2num = 0; // do this like in buildNeighbors
    if(el2->getPartition()==myid){ el2num = ElementMap.at(el2->getId());}
    else                         { el2num = ghostElementMap.at(el2->getId());} // refer to a back column of U
    int id1 = face.getClosureId(0);
    int id2 = face.getClosureId(1);
    const std::vector<int> &cl1 = closures[id1];
    const std::vector<int> &cl2 = closures[id2];
    /*
    printf("id1=%d, id2=%d\n",id1,id2);
    printf("closure::\n");
    for (int m = 0; m < M_s; m++)
      {
	printf("cl1[%d]=%d     cl2[%d]=%d\n",m,cl1[m],m,cl2[m]);
      }
    */
    for(int fc = 0; fc < N_F; fc++){
      for(int j = 0; j < M_s; j++){

	// Map from U to UF
	map[((t*N_F+fc)*2+0)*M_s+j] = (e1*N_F+fc)*N_s+cl1[j];
	map[((t*N_F+fc)*2+1)*M_s+j] = (el2num*N_F+fc)*N_s+cl2[j];

	// Inverse map
	invmap[((e1*N_F+fc)*M_s*N_N+cnt[e1*N_F+fc])*2+0] = cl1[j];
	invmap[((e1*N_F+fc)*M_s*N_N+cnt[e1*N_F+fc])*2+1] = ((t*N_F+fc)*2+0)*M_s+j;

	printf("f=%d,j=%d: just set invmap[[(%d*N_F+%d)*M_s*N_N + %d)*2+0] = %d\n", fc,j,e1,fc,cnt[e1*N_F+fc], cl1[j]);
	printf("f=%d,j=%d: just set invmap[[(%d*N_F+%d)*M_s*N_N + %d)*2+1] = %d\n", fc,j,e1,fc,cnt[e1*N_F+fc], ((t*N_F+fc)*2+0)*M_s+j);
	//PEJ: I think id is always less than N_N for first element (A), geq N_N for second element (B)
	/*
	Alt_FaceFromElem_serial[e1*N_N + id1%N_N] = t;
	printf("Set Alt_FaceFromElem_serial[%d*N_N + %d] = %d\n",e1,id1%N_N, t);
	*/
	if (fc == 0)
	  {
	    //	    Alt_FaceFromElem_serial[e1*N_N + id1%N_N] = t;
	    //printf("Set Alt_FaceFromElem_serial[%d*N_N + %d] = %d\n",e1,id1%N_N, t);
	    Alt_FaceFromElem_serial[e1*(D+N_N) + cnt[e1*N_F+fc]/M_s] = t;
	    if (verbose>0){printf("Set Alt_FaceFromElem_serial[%d*(D+N_N) + %d] = %d\n",e1,cnt[e1*N_F+fc]/M_s, t);}
	  }
	cnt[e1*N_F+fc]++;
	// If this face is not on the boundary (belongs to my
	// partition and physical==0), then you can use the interface
	// to map to an element on my partition	
	if((face.getPhysicalTag()==0)&&(el2->getPartition()==myid)){
	  int e2 = ElementMap.at(el2->getId());
	  invmap[((e2*N_F+fc)*M_s*N_N+cnt[e2*N_F+fc])*2+0] = cl2[j];
	  invmap[((e2*N_F+fc)*M_s*N_N+cnt[e2*N_F+fc])*2+1] = ((t*N_F+fc)*2+1)*M_s+j;
	  //printf("f=%d,j=%d: just set invmap[[(%d*N_F+%d)*M_s*N_N + %d)*2+0] = %d\n", fc,j,e2,fc,cnt[e2*N_F+fc], cl2[j]);
	  //printf("f=%d,j=%d: just set invmap[[(%d*N_F+%d)*M_s*N_N + %d)*2+1] = %d\n", fc,j,e2,fc,cnt[e2*N_F+fc], ((t*N_F+fc)*2+1)*M_s+j);
	  if (fc == 0)
	    {
	      Alt_FaceFromElem_serial[e2*(D+N_N) + cnt[e2*N_F+fc]/M_s] = t;
	      if (verbose>0){printf("Set Alt_FaceFromElem_serial[%d*(D+N_N) + %d] = %d\n",e2,cnt[e2*N_F+fc]/M_s, t);}
	    }
	  cnt[e2*N_F+fc]++;
	  /*
	    Alt_FaceFromElem_serial[e2*N_N + id2%N_N] = t;
	  printf("Set Alt_FaceFromElem_serial[%d*N_N + %d] = %d\n",e2,id2%N_N, t);
	  */
	}
      }
    }
  }
  delete[] cnt;
}
