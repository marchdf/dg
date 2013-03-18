#include "simpleMesh.h"
#include "GmshDefines.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <algorithm>
#include "polynomialBasis.h"
#include "quadratures/Gauss.h"

void simpleMesh::load (const char *fileName)
{
  std::ifstream input (fileName);
  std::string line;
  getline (input, line);
  getline (input, line);
  getline (input, line);
  getline (input, line);
  if (line!="$Nodes")
    printf("invalid file format\n");
  int nnodes;
  input>>nnodes;
  _nodes.resize (3, nnodes);
  for (int i = 0; i < nnodes; i++) {
    int nodeId;
    if (!(input >> nodeId >> _nodes (0, i) >> _nodes (1, i) >> _nodes (2, i)))
      printf("invalid file format\n");
  }
  getline (input, line);
  getline (input, line);
  if (line!="$EndNodes")
    printf("invalid file format\n");
  getline (input, line);
  if (line!="$Elements")
    printf("invalid file format\n");
  int nelements;
  input >> nelements;
  std::vector<int> enodes;
  _elements.resize(MSH_NUM_TYPE);
  getline (input, line);
  for (int i = 0; i < nelements; i++) {
    enodes.resize(0);
    int elementId, elementType, ntags, ptag, num;
    getline (input, line);
    std::istringstream sline (line);
    sline >> elementId >> elementType >> ntags; 
    for (int j = 0; j < ntags; j++) {
      int tag;
      sline >> ((j==1) ? tag : ptag);;
    }
    int idNode;
    while (sline >> idNode) {
      enodes.push_back(idNode-1);
    }
    _elements[elementType].push_back (simpleElement(elementId, ptag, enodes));
  }
  getline (input, line);
  if (line!="$EndElements")
    printf("invalid file format\n");
}

void simpleMesh::buildNormals (int typeInterface, int typeElement, int D)
{
  
  //
  // Initialize the elements, interfaces, nodal functions and closures
  //
  const std::vector<simpleElement> &elements = getElements(typeElement);
  const std::vector<simpleInterface> &interfaces = getInterfaces();
  const polynomialBasis &basis  = *polynomialBases::find (typeElement);  // for the element
  const std::vector<std::vector<int> > &closures = basis.closures;

  // Define some numbers for clarity
  int N_s = elements[0].getNbNodes(); // number of nodes on an element          (i index)
  int N_E = elements.size();          // number of elements                     (e index)
  int N_I = _interfaces.size();       // number of elements                     (i index)
  int N_T = basis.numFaces;           // number of faces per element            (t index)

  //
  // Fonctions de formes et derivees evaluee aux noeuds du maillage
  //
  fullMatrix<double> points2 = basis.points;
  fullMatrix<double> points(points2.size1(),3);
  for(int i = 0; i < points2.size1(); i++){
    for(int j = 0; j < points2.size2(); j++){
      points(i,j) = points2(i,j);
    }
  }
  
  fullMatrix<double> phi (N_s,N_s);    // [phi_j(x_i)]
  fullMatrix<double> dphi(N_s*D,N_s);
  basis.f (points, phi);
  double grads[N_s][3];  
  for(int g = 0; g<N_s; g++){
    basis.df(points(g,0),points(g,1),points(g,2),grads);
    for(int alpha = 0; alpha < D; alpha ++){
      for(int i = 0; i < N_s; i++){
   	dphi(g*D+alpha,i) = grads[i][alpha];  // see paper for indexing p.6
      }	  
    }    
  }

  //
  // Calculate the normals
  //
  _normals.resize(D,N_I);
  _normals.setAll(0.0);
  double* n = new double[D];
  
  // Loop on all the interfaces
  for(int i = 0; i < N_I; i++){

    // Get the interface
    const simpleInterface &face = interfaces[i];

    // Get some information: the element, the closure
    const simpleElement *el = face.getElement(0);
    int clId = face.getClosureId(0);
    const std::vector<int> &cl = closures[clId];

    // Set n to zero
    for(int alpha = 0; alpha < D; alpha++){ n[alpha] = 0;}
    
    // get XYZ coordinates of the element
    fullMatrix<double> XYZNodes (N_s, D);
    for (int j = 0; j < N_s; j++) {
      int inode = el->getNode(j);
      for(int alpha = 0; alpha < D; alpha++){
    	XYZNodes(j, alpha) = _nodes(alpha, inode);
      }
    }
	
    // Jacobian matrix of the element
    fullMatrix<double> Jac(D,N_s*D);
    Jac.gemm(XYZNodes.transpose(),dphi.transpose());

    // Inverse Jacobian matrix 
    fullMatrix<double> invJac(D,N_s*D);
    for(int k = 0; k < N_s; k++ ){
#ifdef ONED
      invJac(0,k) = 1.0/Jac(0,k);
#elif TWOD
      double idet = 1.0/(Jac(0,k*D+0)*Jac(1,k*D+1)-Jac(1,k*D+0)*Jac(0,k*D+1));
      invJac(0,k*D+0) = idet*Jac(1,k*D+1);
      invJac(1,k*D+0) = -1.0*idet*Jac(1,k*D+0);
      invJac(0,k*D+1) = -1.0*idet*Jac(0,k*D+1);
      invJac(1,k*D+1) = idet*Jac(0,k*D+0);
#endif
    }

    // Loop on the nodes of the interface and calculate the normal
    for(int k = 0; k < cl.size(); k++){
      for(int alpha = 0; alpha < D; alpha++){
	for(int a = 0; a < D; a++){ // loop for matrix-matrix mult of invJac and dphi
	  n[alpha] += invJac(a,cl[k]*D+alpha)*dphi(cl[k]*D+a,cl[k]);
	}
      }
    }

    // Normalize and store
    scalar norm2 = 0.0;
    for(int alpha=0; alpha < D; alpha++){ norm2 += n[alpha]*n[alpha];}
    for(int alpha=0; alpha < D; alpha++){
      _normals(alpha,i) = (scalar)n[alpha]/sqrt(norm2);
    }
  } // end loop on interfaces
  delete[] n;
}


void simpleMesh::writeSolution (const fullMatrix<scalar> &solution, int type, const char *filename, const char *name, int step, double time, int append) const
{
  const std::vector<simpleElement> &list = _elements[type];
  if (list.size() != solution.size2())
    printf("bad solution for this element\n");
  std::ofstream output;
  if (append == 1)  output.open(filename,std::ios_base::app);
  else output.open(filename);
  output.precision(20);
  output << "$MeshFormat\n2.1 0 8\n$EndMeshFormat\n";
  output << "$ElementNodeData\n";
  output << "1\n\"" << name << "\"\n";
  output << "1\n" << time << "\n";
  output << "3\n" << step<< "\n1\n" << list.size() << "\n";
  for (int i = 0; i < list.size(); i++) {
    const simpleElement &element = list[i];
    output << element.getId() << " " << element.getNbNodes() <<" ";
    if (element.getNbNodes() != solution.size1())
      printf("bad solution for this element\n");
    for (int j = 0; j < element.getNbNodes(); j++) {
      output << solution (j,i)<< " ";
    }
    output << "\n";
  }
  output << "$EndElementNodeData\n";
}


void simpleMesh::buildLineBoundary(int boundaryType){
  /* Boundary maker in 1D
     boundaryType:
        0 : periodic everywhere
	1 : farfield everywhere
	2 : rflctive everywhere
  */

  // number of boundaries  
  _N_B = 2;
  
  // For the line we build with our mesher, the ends are the first two
  // interfaces
  _boundary = new int[2*_N_B]; // 2 interfaces
  _boundaryIdx = new int[2];
  if(boundaryType==0){ // periodic
    _boundary[0*2+0] = 0;
    _boundary[0*2+1] = 1;
    _boundary[1*2+0] = 1;
    _boundary[1*2+1] = 0;
    _boundaryIdx[0]  = 4; _boundaryIdx[1]  = 4;
  }
  else if(boundaryType==1){ // farfield
    _boundary[0*2+0] = 0;
    _boundary[0*2+1] = 0;
    _boundary[1*2+0] = 1;
    _boundary[1*2+1] = 1;
    _boundaryIdx[0]  = 0;
    _boundaryIdx[0]  = 0; _boundaryIdx[1]  = 4;
  }
  else if(boundaryType==2){ // reflective
    _boundary[0*2+0] = 0;
    _boundary[0*2+1] = 0;
    _boundary[1*2+0] = 1;
    _boundary[1*2+1] = 1;
    _boundaryIdx[0]  = 0; _boundaryIdx[1]  = 0;
  }
}


void simpleMesh::buildSquareBoundary(int M_s, const fullMatrix<scalar> &XYZNodesF, const int D){
  /* Match the interfaces with their partner

     Returns boundaryMap and boundaryidx.
     boundaryMap holds:
     [ [[t1,t2,nodepair][..] [[t1,t2,nodepair][..] [[t1,t2,nodepair][..]]
       periodic              farfield              reflected
     boundaryidx holds:
     the idx of where each of those start

     The boundary type is deduced from the physical of the interface
     (defined when you build the mesh: (1=periodic, 2=farfield,
     3=rflctive)).

     The matching interface pairs are deduced from the (x,y) start and
     end nodes and their normals (opposite and parallel normals means
     the interfaces face each other).
     
  */

  int N_I =_interfaces.size(); // number of interfaces
  
  //
  // Sort each interface according to type
  //
  std::vector<simpleInterface> periodic;
  std::vector<simpleInterface> farfield;
  std::vector<simpleInterface> rflctive;
  std::vector<int> periodicFaceNumber; // these hold the index to the face in the general matrix
  std::vector<int> farfieldFaceNumber;
  std::vector<int> rflctiveFaceNumber;

  for(int i = 0; i < N_I; i++){
    const simpleInterface &face = _interfaces[i];
    int physical = face.getPhysicalTag();
    if     (1==physical){ periodic.push_back(face); periodicFaceNumber.push_back(i); }
    else if(2==physical){ farfield.push_back(face); farfieldFaceNumber.push_back(i); }
    else if(3==physical){ rflctive.push_back(face); rflctiveFaceNumber.push_back(i); }
  }

  //
  // Number of boundary faces
  //
  _N_B = periodic.size() + farfield.size() + rflctive.size();
  
  //
  // Match the boundaries
  // 
  _boundary = new int[2*_N_B];
  _boundaryIdx = new int[2];
  int t = 0;

  // Periodic BC
  for(int t1=0; t1<periodic.size(); t1++){
    // get the face
    const simpleInterface &face1 = periodic[t1];
    int idx1 = periodicFaceNumber[t1];
    // get its end coordinates
    double x0 = MIN(XYZNodesF(0,(idx1*2+0)*D+0),XYZNodesF(1,(idx1*2+0)*D+0));
    double x1 = MAX(XYZNodesF(0,(idx1*2+0)*D+0),XYZNodesF(1,(idx1*2+0)*D+0));
    double y0 = MIN(XYZNodesF(0,(idx1*2+0)*D+1),XYZNodesF(1,(idx1*2+0)*D+1));
    double y1 = MAX(XYZNodesF(0,(idx1*2+0)*D+1),XYZNodesF(1,(idx1*2+0)*D+1));
    double nx1 = _normals(0,idx1);
    double ny1 = _normals(1,idx1);

    for(int t2=0; t2<periodic.size(); t2++){
      if(t2!=t1){ // make sure you don't find yourself
	// get the face
	const simpleInterface &face2 = periodic[t2];
	int idx2 = periodicFaceNumber[t2];
	// get its end coordinates
	double xx0 = MIN(XYZNodesF(0,(idx2*2+0)*D+0),XYZNodesF(1,(idx2*2+0)*D+0));
	double xx1 = MAX(XYZNodesF(0,(idx2*2+0)*D+0),XYZNodesF(1,(idx2*2+0)*D+0));
	double yy0 = MIN(XYZNodesF(0,(idx2*2+0)*D+1),XYZNodesF(1,(idx2*2+0)*D+1));
	double yy1 = MAX(XYZNodesF(0,(idx2*2+0)*D+1),XYZNodesF(1,(idx2*2+0)*D+1));
	double nx2 = _normals(0,idx2);
	double ny2 = _normals(1,idx2);
	
	double eps = 1e-8;
	// Make sure their normals are parallel but opposite (which
	// means the interfaces face each other)
	if(nx1*nx2+ny1*ny2 < 0){ // dot(n1,n2) = |n1|*|n2|*cos(theta) = -1 when // and opposite.
	  // They either have the same x start-end nodes or y start-end nodes
	  if (((fabs(x0-xx0)<eps)&&(fabs(x1-xx1)<eps))||((fabs(y0-yy0)<eps)&&(fabs(y1-yy1)<eps))){
	    _boundary[2*t+0] = idx1;
	    _boundary[2*t+1] = idx2;
	    t++;
	    break; // break out of the loop
	  }
	}
      }
    }
  }
  _boundaryIdx[0] = 2*t;
  
  // Farfield BC
  for(int t1=0; t1<farfield.size(); t1++){
    // get the face
    const simpleInterface &face1 = farfield[t1];
    int idx1 = farfieldFaceNumber[t1];
    _boundary[2*t+0] = idx1;
    _boundary[2*t+1] = idx1;
    t++;
  }
  _boundaryIdx[1] = 2*t;
  
  // Reflective BC
  for(int t1=0; t1<rflctive.size(); t1++){
    // get the face
    const simpleInterface &face1 = rflctive[t1];
    int idx1 = rflctiveFaceNumber[t1];
    _boundary[2*t+0] = idx1;
    _boundary[2*t+1] = idx1;
    t++;
  }

  for(int t=0; t<_N_B; t++){
    printf("B1=%i and B2=%i\n",_boundary[2*t+0], _boundary[2*t+1]);
  }

  // // Match the nodes explicitely
  // for(int b=0; b<_N_B; b++){
  //   int t1 = _boundary[b*2+0];
  //   int t2 = _boundary[b*2+1];
  //   for(int j1 =0; j1<M_s; j1++){
  //     double x1 = XYZNodesF(j1,(t1*2+1)*D+0);
  //     double y1 = XYZNodesF(j1,(t1*2+1)*D+1);
  //     for(int j2 =0; j2<M_s; j2++){
  // 	double x2 = XYZNodesF(j2,(t2*2+0)*D+0);
  // 	double y2 = XYZNodesF(j2,(t2*2+0)*D+1);
  // 	printf("(x1,y1)=(%f,%f) and (x2,y2)=(%f,%f)\n",x1,y1,x2,y2);
  // 	if((fabs(x1-x2)<1e-6)||(fabs(y1-y2)<1e-6)){
  // 	  printf("     j1=%i and j2=%i\n",j1,j2);
  // 	  break;
  // 	}
  //     }
  //   }
  //   printf("\n");
  // }
}


void simpleMesh::buildBoundaryElementShift(int order, const fullMatrix<scalar> &XYZNodesF, const int D, std::map<int,int> &ElementMap)
{
  // Objective: create BoundaryElemShift
  // elem1         | elem2         | xshift | yshift
  // (tgt element)   his neighbor    shifts to bring elem2 next to elem1
  
  _shifts.resize(_N_B,2+D);


  // Most of this code is a repeat of buildPeriodicSquare()
  const std::vector<simpleInterface> &interfaces = getInterfaces();
  
  int N_I = _interfaces.size();       // number of interfaces                     (i index)
  int M_s = order+1;                  // number of nodes on a face

  // Find all boundary interfaces
  int N_B = 0;                        // number of boundary interfaces (b index)
  for(int i = 0; i < N_I; i++){
    const simpleInterface &face = interfaces[i];
    if(face.getPhysicalTag()==1){
      N_B++;
    }
  }
 
  // Get the center of each interface
  double* listInterfaces = new double[3*N_B];
  int counter = 0;
  double xcenmin = 0.0;
  double xcenmax = 0.0;
  double ycenmin = 0.0;
  double ycenmax = 0.0;
  for(int i = 0; i < N_I; i++){
    const simpleInterface &face = interfaces[i];
    if(face.getPhysicalTag()==1){
      int t = i;
      double xcen=0;
      double ycen=0;
      for(int j = 0; j < M_s; j++){
  	xcen += XYZNodesF(j,(t*2+0)*D+0);
  	ycen += XYZNodesF(j,(t*2+0)*D+1);
      }
      listInterfaces[counter*3+0] = t;
      listInterfaces[counter*3+1] = xcen/M_s;
      listInterfaces[counter*3+2] = ycen/M_s;      
      counter++;
      if(xcenmin > xcen/M_s) xcenmin = xcen/M_s;
      if(xcenmax < xcen/M_s) xcenmax = xcen/M_s;
      if(ycenmin > ycen/M_s) ycenmin = ycen/M_s;
      if(ycenmax < ycen/M_s) ycenmax = ycen/M_s;
      //printf("t:%i %i xcen:%f ycen:%f\n",t,face.getPhysicalTag(),xcen/M_s,ycen/M_s);
    }
  }
  double L = xcenmax-xcenmin;
  double H = ycenmax-ycenmin;

  // Match the faces together
  for(int c1 = 0; c1 < N_B; c1++){
    int t1 = listInterfaces[c1*3+0];
    double xcen1 = listInterfaces[c1*3+1];
    double ycen1 = listInterfaces[c1*3+2];
    for(int c2 = 0; c2 < N_B; c2++){
      int t2 = listInterfaces[c2*3+0];
      double xcen2 = listInterfaces[c2*3+1];
      double ycen2 = listInterfaces[c2*3+2];

      // if the interfaces are pairs of each other
      if (((fabs(xcen1-xcen2)<1E-6)&&(fabs(ycen1-ycen2)+1E-6 >= H)) || ((fabs(ycen1-ycen2)<1E-6)&&(fabs(xcen1-xcen2)+1E-6 >= L))){
	const simpleInterface &face1 = interfaces[t1];
	const simpleInterface &face2 = interfaces[t2];
	const simpleElement *el1 = face1.getElement(0); // get the element of face 1
	const simpleElement *el2 = face2.getElement(0); // get the element of face 2
	int el1num = ElementMap[el1->getId()];         // element idx in order of the U matrix
	int el2num = ElementMap[el2->getId()];         // element idx in order of the U matrix
	_shifts(c1,0) = el1num; // el1: target element 
	_shifts(c1,1) = el2num; // el2: his neighbor
	_shifts(c1,2) = xcen1-xcen2; // xshift to bring el2 next to el1
	_shifts(c1,3) = ycen1-ycen2; // yshift to bring el2 next to el1
	//printf("el1num(%i)=%i, el2num(%i)=%i, xshift=%f, yshift=%f\n",el1->getId(),el1num,el2->getId(),el2num, xcen1-xcen2, ycen1-ycen2);
      }
    }
  }

  // Free some stuff
  delete[] listInterfaces;
}

void simpleMesh::buildNeighbors(int N_N, int N_E, std::map<int,int> &ElementMap)
{
  // Objective: find the neighbors of each element
  // get _neighbors, a N_N x N_E vector:
  // | neighbor1 | ...
  // | neighbor2 | ...

  // Allocate neighbors, set to zero
  _neighbors = new int[N_N*N_E]; for(int k=0; k < N_N*N_E; k++){_neighbors[k]=0;}

  int N_I = _interfaces.size();       // number of interfaces                   (i index)

  // A neighbor counter for each element, set to zero
  int* nn = new int[N_E];  for(int k=0; k<N_E; k++){ nn[k]=0;} 
  
  // Loop through each interface in the mesh
  for(int i = 0; i < N_I; i++){
    const simpleInterface &face = _interfaces[i];  // get the interface
    const simpleElement *el1 = face.getElement(0); // get the element on one side
    const simpleElement *el2 = face.getElement(1); // get the element on the other side
    int el1num = ElementMap[el1->getId()];         // element idx in order of the U matrix

    if(el2!=NULL){
      int el2num = ElementMap[el2->getId()];
      //printf("face %2i: el1 = %2i (or %2i) and el2 = %2i (or %2i)\n", i, el1->getId(), el1num, el2->getId(), el2num);
      _neighbors[el1num*N_N+nn[el1num]] = el2num;
      _neighbors[el2num*N_N+nn[el2num]] = el1num;
      nn[el1num]++; // increment the neighbor counter
      nn[el2num]++; // increment the neighbor counter
    }
    else{
      // Find the current interface in the boundaries
      // and get the other interface related to this interface
      int t2 = 0;
      for(int t = 0; t < _N_B; t++) if (_boundary[t*2+0]==i) t2 = _boundary[t*2+1];
      const simpleInterface &face2 = _interfaces[t2];  // get that other interface
      const simpleElement *el3 = face2.getElement(0);  // get the element on one side
      int el3num = ElementMap[el3->getId()];
      //printf("face %2i: el1 = %2i (or %2i) and el3 = %2i (or %2i)\n", i, el1->getId(), el1num, el3->getId(), el3num);
      _neighbors[el1num*N_N+nn[el1num]] = el3num;
      nn[el1num]++; // increment the neighbor counter
    }
  }

  // Free some stuff
  delete[] nn;
}

void simpleMesh::sortNeighbors(const int N_E, const int N_N, const fullMatrix<scalar> XYZCen){
  
  /* Return the neighbors to the left and right and the neighbors
     below and above each element.

     Sort in L-R-D-U order

     Only for a cartesian mesh
  */

  scalar x, y, xn, yn;
  int neighbor;
  int* LRDU = new int[N_N];
  
  // Loop on each element
  for(int e=0; e<N_E; e++){
    // get its cell center
    x = XYZCen(e,0);
    y = XYZCen(e,1);

    // set to -1 so we can check later which ones were accessed
    for(int nn=0; nn<N_N; nn++) LRDU[nn] = -1; 
    
    // Loop on all the neighbors
    for(int nn=0; nn<N_N; nn++){
      neighbor = _neighbors[e*N_N+nn]; // get a neighbor element
      xn = XYZCen(neighbor,0);        // get the centroid
      yn = XYZCen(neighbor,1);

      // for periodic boundary condition
      // Get potentional coordinate shift for the neighbor if this
      // element is on the boundary
      bool flag = false; int bidx = 0;
      for(int b=0; b<_N_B; b++)
	if ((e==_shifts(b,0))&&(neighbor==_shifts(b,1))){ flag = true; bidx = b; break;}
      if(flag){
	xn = xn+_shifts(bidx,2+0);
	yn = yn+_shifts(bidx,2+1);
      }
      
      // Relation to current element
      if     ((fabs(yn-y)<1e-9)&&(x>xn)){  LRDU[0]=neighbor; /*printf("found L\n");*/}// left
      else if((fabs(yn-y)<1e-9)&&(x<xn)){  LRDU[1]=neighbor; /*printf("found R\n");*/} // right
      else if((fabs(xn-x)<1e-9)&&(y>yn)){  LRDU[2]=neighbor; /*printf("found D\n");*/} // below
      else if((fabs(xn-x)<1e-9)&&(y<yn)){  LRDU[3]=neighbor; /*printf("found U\n");*/} // above
    }//loop on neighbors

    // For the ones that were not allocated, assume it's because we
    // are dealing with a farfield/reflective boundary so we just copy
    // the current element to his neighbor
    for(int nn=0; nn<N_N; nn++) if(LRDU[nn]==-1) LRDU[nn] = e;

    // Sort the neighbors
    for(int nn=0; nn<N_N; nn++) _neighbors[e*N_N+nn] = LRDU[nn];
    
  }// loop on elements

  delete[] LRDU;
}


simpleInterface::simpleInterface(int physicalTag)
{
  _elements[0] = NULL;
  _elements[1] = NULL;
  _physicalTag = physicalTag;
  _closureId[0]= -1;
  _closureId[1]= -1;
}


void simpleInterface::BuildInterfaces(simpleMesh &mesh, std::vector<simpleInterface> &interfaces, int typeInterface, int typeElement, int nsides)
{
  std::map<std::vector<int>, simpleInterface> interfacesMap;
  // pre-existing interfaces
  const std::vector<simpleElement> preElements = mesh.getElements(typeInterface);
  //  printf("preElements size %i\n", preElements.size());
  for (int i = 0; i < preElements.size(); i++) {
    const simpleElement &el = preElements[i];
    std::vector<int> nodes;
    //    printf("element %i\n", i);
    for (int k = 0; k < el.getNbNodes(); k++) {
      nodes.push_back(el.getNode(k));
      //      printf("     node %i\n",el.getNode(k)+1);
    }
    std::sort(nodes.begin(), nodes.end());
    interfacesMap[nodes] = simpleInterface(el.getPhysicalTag());
  }
  //  printf("Map size is %i after the preElement analysis\n", interfacesMap.size());
  // other interfaces
  const std::vector<simpleElement> &elements = mesh.getElements(typeElement);
  const polynomialBasis &basis = *polynomialBases::find(typeElement);
  const std::vector<std::vector<int> > &closures = basis.closures;
  std::vector<int> nodes;
  for (int i = 0; i < elements.size(); i++) {
    const simpleElement &el = elements[i];
    // printf("element %i\n", i);
    // printf("   id of the element is %i\n",(el).getId());
    // printf("   num of nodes is %i\n",(el).getNbNodes());
    // printf("   list of nodes:");
    // for(int p=0; p< el.getNbNodes();p++){
    //   printf(" %i",el.getNode(p));
    // }
    // printf("\n");
    for (int j = 0; j < closures.size(); j++) {
      nodes.clear();
      const std::vector<int> &cl = closures[j];
      for (int k = 0; k < cl.size(); k++) {
        nodes.push_back(el.getNode(cl[k]));
  	//	printf("   node %i\n", nodes[k]);
      }
      std::sort(nodes.begin(), nodes.end());
      simpleInterface &interfaceFound = interfacesMap[nodes];
      if (interfaceFound._elements[0] == NULL) {
        interfaceFound._elements[0] = &el;
  	interfaceFound._elements[1] = NULL;
        interfaceFound._closureId[0] = j;
  	interfaceFound._closureId[1] = -1;
  	//	printf("This has to be 8 times\n");
      }
      else if (interfaceFound._elements[0] != &el) {
        if (interfaceFound._elements[1] == NULL) {
          interfaceFound._elements[1] = &el;
          interfaceFound._closureId[1] = j+nsides; // WHY?! j+3 in 2D (triangle)
  	  //	  printf("This has to be 4 times\n");
        }
  	else if (interfaceFound._elements[1] != &el) {
          printf("error in interfaces creation !!!\n");
          exit(1);
        }
      }
    }
  }
  interfaces.clear();
  interfaces.reserve(interfacesMap.size());
  for (std::map<std::vector<int>, simpleInterface>::iterator it = interfacesMap.begin(); it != interfacesMap.end(); ++it) {
    interfaces.push_back(it->second);
  }
}

// Return a matrix which is N_E x D containing the element centroids
fullMatrix<scalar> simpleMesh::getElementCentroids(const int N_E, const int D, const int ncorners, const fullMatrix<scalar> XYZNodes){

  fullMatrix<scalar> XYZCen(N_E,D);
  scalar cen =0;
  
  for(int e = 0; e < N_E; e++){
    for(int alpha = 0; alpha < D; alpha++){
      cen = 0;
      for(int k = 0; k < ncorners; k++){
	cen += XYZNodes(k,e*D+alpha);	
      }
      XYZCen(e,alpha) = cen/ncorners;
    }    
  }
  return XYZCen;
}

bool simpleMesh::iscartesian(std::string typeElement, const int elem_type){

  // default is unstructured mesh
  bool cartesian = false;
  
  if (typeElement=="qua"){
    const std::vector<simpleElement> &elements = getElements(elem_type);
    const simpleElement &el = elements[0];
    scalar x0 = (scalar)_nodes(0,el.getNode(0));
    scalar y0 = (scalar)_nodes(1,el.getNode(0));
    scalar x1 = (scalar)_nodes(0,el.getNode(1));
    scalar y1 = (scalar)_nodes(1,el.getNode(1));
    scalar x2 = (scalar)_nodes(0,el.getNode(2));
    scalar y2 = (scalar)_nodes(1,el.getNode(2));
    scalar x3 = (scalar)_nodes(0,el.getNode(3));
    scalar y3 = (scalar)_nodes(1,el.getNode(3));

    // Check if the vectors b0 and b1 and b2 are perpendicular
    // b1 = [x1-x0;y1-y0], b1 = [x2-x1;y2-y1] and b2 = [x3-x2;y3-y2]
    if((fabs((x1-x0)*(x2-x1)+(y1-y0)*(y2-y1))<1e-9)&&(fabs((x2-x1)*(x3-x2)+(y2-y1)*(y3-y2))<1e-9)){
      cartesian = true;
    }
  }

  return cartesian;
}


simpleMesh::~simpleMesh(){
  //delete[] _boundary;
  //delete _boundaryIdx;
}
