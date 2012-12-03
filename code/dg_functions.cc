#include <dg_functions.h>
void dg_jac_elements_fast(const int N_G, const int N_E, const int D, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &dphi, fullMatrix<scalar> &J){
  fullMatrix<scalar> Jac(N_E*D,N_G*D);

  Jac.gemm(XYZNodes.transpose(),dphi.transpose());
  scalar det = 0.0;
  
  for(int e = 0; e < N_E; e++){
    det = sqrt(Jac(e,0)*Jac(e,0));    //Jac(e*D+0,0)*Jac(e*D+1,1)-Jac(e*D+1,0)*Jac(e*D+0,1);
    J(e,0) =  fabs(det);
  }
}

void dg_jacobians_elements(const int N_G, const int N_E, const int D, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &dphi, fullMatrix<scalar> &Jac, fullMatrix<scalar> &invJac, fullMatrix<scalar> &J, fullMatrix<scalar> &invJ){

  Jac.gemm(XYZNodes.transpose(),dphi.transpose());
  scalar det = 0.0;
  
  for(int e = 0; e < N_E; e++){
    det = sqrt(Jac(e,0)*Jac(e,0));    //Jac(e*D+0,0)*Jac(e*D+1,1)-Jac(e*D+1,0)*Jac(e*D+0,1);
    J(e,0)    =  fabs(det);
    invJ(e,0) =  fabs(1.0/det);
    for(int g = 0; g < N_G; g++){
      //printf("Jac of e:%i and g:%i = %f\n",e,g,Jac(e,g));
      invJac(e,g) =  1.0/Jac(e,g);
    }
  }
}


void dg_jacobians_face(const int M_T, const int D, fullMatrix<scalar> &XYZNodesF, fullMatrix<scalar> &dpsi, fullMatrix<scalar> &JacF, fullMatrix<scalar> &JF, fullMatrix<scalar> &invJF){

  JacF.gemm(XYZNodesF.transpose(),dpsi.transpose());
  for(int t = 0; t < M_T; t++){
    for(int d = 0; d < 2; d++){
      scalar det = 0.0;
      for(int alpha = 0; alpha < D; alpha ++){ det += sqrt(JacF((t*2+d)*D+alpha,0)* JacF((t*2+d)*D+alpha,0));}
      if(D==1) det = 1; // in 1D set JF to 1
      JF(t*2+d,0)    =  det;
      invJF(t*2+d,0) =  1.0/det;
    }
  }
}


void dg_inverse_mass_matrix(const int order, const int elem_type, const std::string getElemType, const int N_s, const int N_E, const int D, fullMatrix<scalar> &XYZNodes, scalar* Minv){

  const polynomialBasis *basis  = polynomialBases::find (elem_type);  // for the element
  fullMatrix<double> points, weight;
  if(getElemType == "lin") gaussIntegration::getLine(order*2+1, points, weight);

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
  

  // Calculate the jacobians
  fullMatrix<scalar> J(N_E,1);            // determinant of the Jacobian
  dg_jac_elements_fast(N_G, N_E, D, XYZNodes, dphi, J);

  for(int e = 0; e < N_E; e++){
    // Init de psi_w_Jac qui va etre construit pour chaque element
    fullMatrix<scalar> phi_w_Jac = phi_w;
    // Construction de la matrice de masse: M(i,j) = \int psi_i psi_j
    fullMatrix<scalar> M(N_s,N_s);
    scalar jac = J(e,0);
    for(int g = 0; g < N_G; g++){
      for(int i = 0; i < N_s; i++){
        phi_w_Jac(g,i) *= jac;
      }
    }
    M.gemm(phi.transpose(),phi_w_Jac);
    
    // Inverser la matrice de masse.
    fullMatrix<scalar> M_inv(N_s,N_s);
    M.invert(M_inv);
    
    for(int k = 0; k < N_s; k++){
      for(int i = 0; i < N_s; i++){
	Minv[(e*N_s+k)*N_s+i] = M_inv(i,k);
      }
    }
  }
}

void dg_mappings(const int M_s, const int M_T, const int N_F, const int N_s, const int N_E, const std::vector<simpleInterface> &interfaces, std::map<int,int> &ElementMap, const std::vector<std::vector<int> > &closures, int* map, int* invmap){

  for(int k = 0; k < M_s*M_T*N_F*2 ; k++) { map[k] = -1;}
  for(int k = 0; k < N_s*N_E*N_F*2 ; k++) { invmap[k] = -1;}

  for(int t = 0; t < M_T; t++){
    const simpleInterface &face = interfaces[t];    // get the interface
    for(int d = 0; d < 2; d++){
      const simpleElement *el = face.getElement(d); // get the element to the right/left
      if(el!=NULL){                                 // as long as there is a element to the right/left
  	int id = face.getClosureId(d);
  	const std::vector<int> &cl = closures[id];
	for(int fc = 0; fc < N_F; fc++){
	  for(int j = 0; j < M_s; j++){
	    map[((t*N_F+fc)*2+d)*M_s+j] = (ElementMap[el->getId()]*N_F+fc)*N_s+cl[j];
	    if (invmap[((ElementMap[el->getId()]*N_F+fc)*2+0)*N_s+cl[j]] == -1){
	      invmap[((ElementMap[el->getId()]*N_F+fc)*2+0)*N_s+cl[j]] = ((t*N_F+fc)*2+d)*M_s+j;
	    }
	    else if(invmap[((ElementMap[el->getId()]*N_F+fc)*2+0)*N_s+cl[j]] != -1){
	      invmap[((ElementMap[el->getId()]*N_F+fc)*2+1)*N_s+cl[j]] = ((t*N_F+fc)*2+d)*M_s+j;
	    }
	  }
	}
      }
    }
  }    
}


