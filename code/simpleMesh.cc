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
    if (D == 1){ for(int k = 0; k < N_s; k++ ){ invJac(0,k) = 1.0/Jac(0,k);}}
    else if(D==2){
      for(int k = 0; k < N_s; k++ ){
    	double idet = 1.0/(Jac(0,k*D+0)*Jac(1,k*D+1)-Jac(1,k*D+0)*Jac(0,k*D+1));
    	invJac(0,k*D+0) = idet*Jac(1,k*D+1);
    	invJac(1,k*D+0) = -1.0*idet*Jac(1,k*D+0);
    	invJac(0,k*D+1) = -1.0*idet*Jac(0,k*D+1);
    	invJac(1,k*D+1) = idet*Jac(0,k*D+0);
      }
    }

    // // Loop on the nodes of the interface and calculate the normal
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

void simpleMesh::buildPeriodicSquare(int order, const fullMatrix<scalar> &XYZNodesF, const int D)
{
  // Objective: match the interfaces with their periodic partner
  const std::vector<simpleInterface> &interfaces = getInterfaces();
  
  int N_I = _interfaces.size();       // number of interfaces                     (i index)
  int M_s = order+1;                  // number of nodes on a face
  
  // Find all boundary interfaces
  _N_B = 0;                        // number of boundary interfaces (b index)
  for(int i = 0; i < N_I; i++){
    const simpleInterface &face = interfaces[i];
    if(face.getPhysicalTag()==1){
      _N_B++;
    }
  }
  printf("M_B %i\n",_N_B);
  
  // Get the center of each interface
  double* listInterfaces = new double[3*_N_B];
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
  _boundary = new int[2*_N_B];
  for(int c1 = 0; c1 < _N_B; c1++){
    int t1 = listInterfaces[c1*3+0];
    double xcen1 = listInterfaces[c1*3+1];
    double ycen1 = listInterfaces[c1*3+2];
    for(int c2 = 0; c2 < _N_B; c2++){
      int t2 = listInterfaces[c2*3+0];
      double xcen2 = listInterfaces[c2*3+1];
      double ycen2 = listInterfaces[c2*3+2];

      if ((fabs(xcen1-xcen2)<1E-6)&&(fabs(ycen1-ycen2)+1E-6 >= H)){_boundary[c1*2+0] = t1; _boundary[c1*2+1] = t2;}// t1=t2
      if ((fabs(ycen1-ycen2)<1E-6)&&(fabs(xcen1-xcen2)+1E-6 >= L)){_boundary[c1*2+0] = t1; _boundary[c1*2+1] = t2;}// t1=t2
    }
  }

  delete[] listInterfaces;
}

void simpleMesh::buildPeriodicLine()
{
  // Objective: match the interfaces with their periodic partner
  const std::vector<simpleInterface> &interfaces = getInterfaces();

  // For the line we build with our mesher, the ends are the first two
  // interfaces
  _N_B = 2;
  _boundary = new int[2*_N_B]; // 2 interfaces
  _boundary[0*2+0] = 0;
  _boundary[0*2+1] = 1;
  _boundary[1*2+0] = 1;
  _boundary[1*2+1] = 0;
  
}

void simpleMesh::buildFarfield()
{
  // Objective: match the interfaces with their farfield partner
  const std::vector<simpleInterface> &interfaces = getInterfaces();
  
  int N_I = _interfaces.size();       // number of interfaces                     (i index)

  // Find all boundary interfaces
  _N_B = 0;                        // number of boundary interfaces (b index)
  for(int i = 0; i < N_I; i++){
    const simpleInterface &face = interfaces[i];
    if(face.getElement(1)==NULL){ // if the inferface only has one element
      _N_B++;
    }
  }
  printf("M_B %i\n",_N_B);
  
  // Get the center of each interface
  double* listInterfaces = new double[_N_B];
  int counter = 0;
  for(int i = 0; i < N_I; i++){
    const simpleInterface &face = interfaces[i];
    if(face.getElement(1)==NULL){ // if the inferface only has one element
      int t = i;
      listInterfaces[counter] = t;
      counter++;
    }
  }

  // Match the faces together
  _boundary = new int[2*_N_B];
  for(int t = 0; t < _N_B; t++){
    _boundary[t*2+0] = listInterfaces[t];
    _boundary[t*2+1] = listInterfaces[t];
  }

  delete[] listInterfaces;
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

