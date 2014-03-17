#include "simpleMesh.h"
#include "GmshDefines.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <algorithm>
#include "polynomialBasis.h"
#include "quadratures/Gauss.h"
#include <stdexcept>      // std::out_of_range

bool pairPeriodic(const fullMatrix<double> &meshNodes, std::vector<int> nodes1, std::vector<int> nodes2){

  /* Given a list of nodes on an interface and their coordinates,
     figure out if they are periodic pairs.

     Return true if periodic pairs.

     Used by BuildInterfaces. */
     
  bool paired = false;
  double eps = 1e-10;
  
#ifdef ONED
  // Get interface x position of interfaces
  scalar x1 = meshNodes(0,nodes1[0]);
  scalar x2 = meshNodes(0,nodes2[0]);

  // Make sure not to find myself
  if      (fabs(x1-x2)<eps){return paired;}
  // if we are distant (and our physicals==3) then we are paired
  else if (fabs(x1-x2)>eps){paired = true;}

#elif TWOD

  // Get start and end nodes of the interfaces
  scalar xs1 = meshNodes(0,nodes1[0]);
  scalar xe1 = meshNodes(0,nodes1[0]);
  scalar ys1 = meshNodes(1,nodes1[0]);
  scalar ye1 = meshNodes(1,nodes1[0]);
  scalar xs2 = meshNodes(0,nodes2[0]);
  scalar xe2 = meshNodes(0,nodes2[0]);
  scalar ys2 = meshNodes(1,nodes2[0]);
  scalar ye2 = meshNodes(1,nodes2[0]);
  for(int j=1; j<nodes1.size(); j++){
    xs1 = MIN(xs1,meshNodes(0,nodes1[j]));
    xe1 = MAX(xe1,meshNodes(0,nodes1[j]));
    ys1 = MIN(ys1,meshNodes(1,nodes1[j]));
    ye1 = MAX(ye1,meshNodes(1,nodes1[j]));
    xs2 = MIN(xs2,meshNodes(0,nodes2[j]));
    xe2 = MAX(xe2,meshNodes(0,nodes2[j]));
    ys2 = MIN(ys2,meshNodes(1,nodes2[j]));
    ye2 = MAX(ye2,meshNodes(1,nodes2[j]));
  }

  // Make sure not to find myself
  if ((fabs(xs1-xs2)<eps)&&(fabs(xe1-xe2)<eps)&&(fabs(ys1-ys2)<eps)&&(fabs(ye1-ye2)<eps)){ return paired;}
  // vertical interface
  else if (fabs(ys1-ye1)>eps){
    // if their start and end nodes match, pair them
    if ((fabs(ys1-ys2)<eps)&&(fabs(ye1-ye2)<eps)){ paired = true;}
  }
  // horizontal interface
  else if (fabs(xs1-xe1)>eps){
    // if their start and end nodes match, pair them
    if ((fabs(xs1-xs2)<eps)&&(fabs(xe1-xe2)<eps)){ paired = true;}
  }	  
#endif
  
  return paired;  
}


int uniqueTag(int id1, int id2, int N=0){
  /* Given two numbers (the global id of el1 and el2 of an interface),
     create a unique number. Ideally I would use a perfect hash. But
     for that I would need a 64-bit integer to represent all pairs of
     32-bit integers. For a good discussion on this see
     http://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way
     http://stackoverflow.com/questions/682438/hash-function-providing-unique-uint-from-an-integer-coordinate-pair
     http://stackoverflow.com/questions/11786635/fast-bi-directional-hash-of-two-integers-in-c */
  
  int k1 = MIN(id1,id2); 
  int k2 = MAX(id1,id2);

  //int tag = cantor(k1,k2); // this will overflow

  // I decided to go with this solutions. N is the total number of
  // pairs I have. So the condition in
  // http://stackoverflow.com/questions/682438/hash-function-providing-unique-uint-from-an-integer-coordinate-pair
  // is no longer a restriction. Collisions might happen but never for
  // the same processors communicating.
  int tag = N*k1+k2; // This will most likely stay in the 32-bit domain
  
  return tag;
}

int cantor(int k1, int k2){
  /* Cantor pairing function. This will overflow when for
     (k1,k2)=(65535, 65535). You can represent all 16-bit integers
     using 32-bit integers but you can't go further than that without
     using a 64-bit integer.*/
  return 0.5*(k1+k2)*(k1+k2+1)+k2;
}
    
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
  _elements.resize(MSH_NUM_TYPE);      // holds elements in my partition
  _otherElements.resize(MSH_NUM_TYPE); // holds elements in other partitions
  getline (input, line);
  for (int i = 0; i < nelements; i++) {
    enodes.resize(0);
    int elementId, elementType, ntags, ptag, num, partition=1;
    getline (input, line);
    std::istringstream sline (line);
    sline >> elementId >> elementType >> ntags; 
    for (int j = 0; j < ntags; j++) {
      int tag;
      if      (j==0) sline >> ptag;      // physical tag of element
      else if (j==3) sline >> partition; // main partition of element
      else           sline >> tag;
      //printf("Id=%i: j=%i, tag=%i, ptag=%i, partition=%i \n",elementId, j, tag,ptag, partition);
    }
    int idNode;
    while (sline >> idNode) {
      enodes.push_back(idNode-1);
    }
    // Exit if the partition number is larger than the total number of processors
    if(_numprocs < partition){printf("Not enough processors for the mesh partitions. Exiting\n"); exit(1);}
    // Only store the element if it's in my partition
    if(_myid == partition-1) _elements[elementType].push_back (simpleElement(elementId, ptag, partition-1, enodes));
    // Otherwise store the element in the otherElements
    else _otherElements[elementType].push_back (simpleElement(elementId, ptag, partition-1, enodes));
  }
  getline (input, line);
  if (line!="$EndElements")
    printf("invalid file format\n");
}

void simpleMesh::buildNormals (int typeInterface, int typeElement)
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


void simpleMesh::writeSolution (const fullMatrix<scalar> &solution, int type, std::string filename, std::string name, int step, double time, int append) const
{

#ifdef USE_MPI
  char numstr[21]; // enough to hold all numbers up to 64-bits
  sprintf(numstr, "%d", _myid);
  filename += numstr;
#endif
  
  const std::vector<simpleElement> &list = _elements[type];
  if (list.size() != solution.size2())
    printf("bad solution for this element\n");
  std::ofstream output;
  if (append == 1)  output.open(filename.c_str(),std::ios_base::app);
  else output.open(filename.c_str());
  output.precision(20);
  output << "$MeshFormat\n2.1 0 8\n$EndMeshFormat\n";
  output << "$ElementNodeData\n";
  output << "1\n\"" << name << "\"\n";
  output << "1\n" << time << "\n";
  output << "4\n" << step<< "\n1\n" << list.size() << "\n" << _myid << "\n";
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

void simpleMesh::buildBoundary(){
  /*
    This is just going to get us a list of the interfaces where the
    boundaries are not farfield or periodic. For example, it will have
    the list of reflective BC (and other, more complicated, ones if we
    want).
    
    Returns boundaryMap and boundaryidx.
    boundaryMap holds:
    [ [t1 ..]     [t1 ..]]
      reflective  otherone
    boundaryidx holds:
    the idx of where each of those start
    
    The boundary type is deduced from the physical of the interface
    (defined when you build the mesh: (3=rflctive)).
  */

  int N_I =_interfaces.size(); // number of interfaces
  
  //
  // Sort each interface according to type
  //
  std::vector<simpleInterface> rflctive;
  std::vector<simpleInterface> otherone;

  std::vector<int> rflctiveFaceNumber; // these hold the index to the face in the general matrix
  std::vector<int> otheroneFaceNumber; 

  for(int i = 0; i < N_I; i++){
    const simpleInterface &face = _interfaces[i];
    int physical = face.getPhysicalTag();
    if     (3==physical){ rflctive.push_back(face); rflctiveFaceNumber.push_back(i);}
    else if(5==physical){ otherone.push_back(face); otheroneFaceNumber.push_back(i);}
  }

  //
  // Number of boundary faces
  //
  _N_B = rflctive.size() + otherone.size();

  //
  // Match the boundaries
  // 
  _boundary = new int[_N_B];
  _boundaryIdx = new int[1];
  int t = 0;

  // Reflective BC
  for(int t1=0; t1<rflctive.size(); t1++){
    // get the face
    const simpleInterface &face1 = rflctive[t1];
    int idx1 = rflctiveFaceNumber[t1];
    _boundary[t] = idx1;
    t++;
  }
  _boundaryIdx[0] = t;
  
  // Other BC
  for(int t1=0; t1<otherone.size(); t1++){
    // get the face
    const simpleInterface &face1 = otherone[t1];
    int idx1 = otheroneFaceNumber[t1];
    _boundary[t] = idx1;
    t++;
  }
  _boundaryIdx[0] = t;
}


void simpleMesh::buildBoundaryElementShift1D(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes){ 
  // Objective: create BoundaryElemShift
  // elem1         | elem2         | xshift 
  // (tgt element)   his neighbor    shifts to bring elem2 next to elem1

  // In a 1D straight line, the leftmost element is 0 and the rightmost element is N_E-1
  // The shift to bring them next to each other is equal to the length of the domain.
  scalar xmin = XYZNodes(0,0);
  scalar xmax = XYZNodes(0,0);
  scalar x = 0;
  for(int e=0; e<N_E; e++){
    for(int i=0; i<N_s; i++){
      x = XYZNodes(i,e);
      if(xmin>x) xmin = x;
      if(xmax<x) xmax = x;
    }
  }
  scalar L = xmax - xmin;

  // Build shifts
  _shifts.resize(_N_B,2+D);
  _shifts(0,0) = 0;
  _shifts(0,1) = N_E-1;
  _shifts(0,2) =-L; // xshift to bring el2 next to el1
  _shifts(1,0) = N_E-1;
  _shifts(1,1) = 0;
  _shifts(1,2) = L;

}
  
void simpleMesh::buildBoundaryElementShift2D(int order, const fullMatrix<scalar> &XYZNodesF, std::map<int,int> &ElementMap)
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

void simpleMesh::buildNeighbors(int N_N, int N_E)
{
  // Objective: find the neighbors of each element
  // get _neighbors, a N_N x N_E vector:
  // | neighbor1 | ...
  // | neighbor2 | ...

  // Allocate neighbors, set to zero
  _neighbors = new int[N_N*N_E]; for(int k=0; k < N_N*N_E; k++){_neighbors[k]=0;}

  int N_I = _interfaces.size();       // number of interfaces                   (i index)
  double eps = 1e-9; // used for comparisons of smallness
  
  // A neighbor counter for each element, set to zero
  int* nn = new int[N_E];  for(int k=0; k<N_E; k++){ nn[k]=0;} 

  // Loop through each interface in the mesh
  for(int i = 0; i < N_I; i++){
    const simpleInterface &face = _interfaces[i];  // get the interface
    const simpleElement *el1 = face.getElement(0); // get the element on one side
    const simpleElement *el2 = face.getElement(1); // get the element on the other side
    int el1num = _elementMap.at(el1->getId());     // element idx in order of the U matrix
    int el2num=0; 
    if(el2->getPartition()==_myid){el2num = _elementMap.at(el2->getId());}
    else                          {el2num = _ghostElementMap.at(el2->getId());}


    // If el1num = el2num, we are at some kind of boundary
    if (el2num == el1num){
      int physical = face.getPhysicalTag();
      // give elnum2 the negative value of the physical. The idea is
      // that when using the neighbors in the limiting procedure, we can
      // know which boundary we are at and use the appropriate nodal
      // values for the reconstruction
      if ((physical == 2)||(physical == 3)||(physical == 4)){el2num = -physical;}
    }
    
    // Not a cartesian mesh
    if(!_cartesian){
      _neighbors[el1num*N_N+nn[el1num]] = el2num; nn[el1num]++;
      // do the reverse if el2 belongs to me
      if (el2num>=0){ if(el2->getPartition()==_myid){ _neighbors[el2num*N_N+nn[el2num]] = el1num; nn[el2num]++;}}
    }
    else if (_cartesian){ // sort in LRDU order
      double nx = _normals(0,i);
#ifdef ONED
      double ny = 0;
#elif TWOD
      double ny = _normals(1,i);
#endif
      // printf("e1=%i, e2=%i, nx=%e, ny=%e\n",el1->getId(),el2->getId(),nx,ny);
      // normal pointing to the left: el2 is at left of el1
      if     ((fabs(ny)<eps)&&(nx<0)){
	_neighbors[el1num*N_N+0] = el2num;
	// do the reverse if el2 belongs to me (and it's not el1num, ie el2num >=0 )
	if (el2num>=0){ if(el2->getPartition()==_myid){_neighbors[el2num*N_N+1] = el1num;}}
      }

      // normal pointing to the right: el2 is at right of el1
      else if((fabs(ny)<eps)&&(nx>0)){
	_neighbors[el1num*N_N+1] = el2num;
	if (el2num>=0){ if(el2->getPartition()==_myid){ _neighbors[el2num*N_N+0] = el1num;}}
      }

      // normal pointing down: el2 is below el1
      else if((fabs(nx)<eps)&&(ny<0)){
	_neighbors[el1num*N_N+2] = el2num;
	if (el2num>=0){ if(el2->getPartition()==_myid){ _neighbors[el2num*N_N+3] = el1num;}}
      }

      // normal pointing up: el2 is above el1
      else if((fabs(nx)<eps)&&(ny>0)){
	_neighbors[el1num*N_N+3] = el2num;
	if (el2num>=0){ if(el2->getPartition()==_myid){ _neighbors[el2num*N_N+2] = el1num;}}
      }
      else{
	printf("Error in neighbor creation!!\n");
	exit(1);
      }
    } // end if cartesian
  } // end loop on interfaces

  // Free some stuff
  delete[] nn;
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
  /* This function builds all the interfaces. For each interface, it
     will find both elements on the right and left. If it's on a
     boundary, it will use the BC information to find the correct
     neighboring element.

     For MPI processes, there a bunch of conditions to take care of so
     things might get a bit convoluted. Hang on to your hats.*/

  std::map<std::vector<int>, simpleInterface> interfacesMap;
  //
  // Read the pre-existing interfaces in the mesh (my partition)
  //
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

  //
  // Now loop through the elements in the mesh and find all the other
  // interfaces (in my partition)
  //
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
      
      // If found in the map, use it. If not, add it to the map
      simpleInterface &interfaceFound = interfacesMap[nodes];
      if (interfaceFound._elements[0] == NULL) {
        interfaceFound._elements[0] = &el;
        interfaceFound._closureId[0] = j;
	//farfield or reflective BC, copy my element to my neighbor
	if ((interfaceFound.getPhysicalTag()==2)||(interfaceFound.getPhysicalTag()==3)||(interfaceFound.getPhysicalTag()==4)){
	  interfaceFound._elements[1] = &el;
	  interfaceFound._closureId[1] = j;
	}
	else{	  
	  interfaceFound._elements[1] = NULL;
	  interfaceFound._closureId[1] = -1;
	}
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

  //
  // Now account for the periodic BC which are on my partition
  //
  const fullMatrix<double> &meshNodes = mesh.getNodes();  
  for (std::map<std::vector<int>, simpleInterface>::iterator it = interfacesMap.begin(); it != interfacesMap.end(); ++it) {
    std::vector<int> nodes1 = it -> first;
    simpleInterface & interface1 = it -> second;

    // If it's a periodic BC and we haven't found a partner element yet
    if ((interface1.getPhysicalTag()==1)&&(interface1.getElement(1)==NULL)){
  
      // Now loop on all the other interfaces to find the match
      for (std::map<std::vector<int>, simpleInterface>::iterator it2 = interfacesMap.begin(); it2 != interfacesMap.end(); ++it2) {
  	std::vector<int> nodes2 = it2 -> first;
  	simpleInterface & interface2 = it2 -> second;

  	// If it's a periodic BC and we haven't found a partner element yet
  	if ((interface2.getPhysicalTag()==1)&&(interface2.getElement(1)==NULL)){

  	  // Check if they are paired
  	  bool paired = pairPeriodic(meshNodes, nodes1, nodes2);
  	  if(paired){
  	    interface1._elements[1] = interface2.getElement(0);
  	    interface1._closureId[1] = interface2.getClosureId(0)+nsides;
  	    interface2._elements[1] = interface1.getElement(0);
  	    interface2._closureId[1] = interface1.getClosureId(0)+nsides;
  	  }
  	} // end second if condition on periodic
      } // end second loop on interfaceMap
    } // end first if condition on periodic
  }// end first loop on interfaceMap


  //
  // Classify all periodic BC that are not in my partition
  //
  std::map<std::vector<int>, simpleInterface> ghostPeriodicMap;
  const std::vector<simpleElement> ghostPreElements = mesh.getOtherElements(typeInterface);
  for(int i = 0; i < ghostPreElements.size();i++){
    const simpleElement &el = ghostPreElements[i];
    if(el.getPhysicalTag()==1){ // periodic boundary
      nodes.clear();
      for (int k = 0; k < el.getNbNodes(); k++) {
  	nodes.push_back(el.getNode(k));
      }
      std::sort(nodes.begin(), nodes.end());
      ghostPeriodicMap[nodes] = simpleInterface(el.getPhysicalTag());
    } // end if periodic
  } // end loop on other interfaces

  //
  // Now take care of all the inner boundaries (the interfaces
  // directly connected to elements on another partition). Very
  // similar to first interface-element matching. Also deal with all
  // the periodic boundaries which have pairs in other partitions
  //
  const std::vector<simpleElement> &otherElements = mesh.getOtherElements(typeElement);
  for (int i = 0; i < otherElements.size(); i++) {
    const simpleElement &el = otherElements[i];
    for (int j = 0; j < closures.size(); j++) {
      nodes.clear();
      const std::vector<int> &cl = closures[j];
      for (int k = 0; k < cl.size(); k++) {
        nodes.push_back(el.getNode(cl[k]));
      }
      std::sort(nodes.begin(), nodes.end());
      
      // If found in the map, use it. If not, move on (do not add to the map)
      try{
  	simpleInterface &interfaceFound = interfacesMap.at(nodes);
  	// If this interface has a first element but not a second element
  	if ((interfaceFound._elements[0] != NULL)&&(interfaceFound._elements[1] == NULL)) {
  	  interfaceFound._elements[1] = &el;
  	  interfaceFound._closureId[1] = j+nsides;
  	}
      }
      catch(const std::out_of_range &oor){
  	// See if this interface is in the ghostPeriodicMap
  	try{
  	  simpleInterface &periodicFound = ghostPeriodicMap.at(nodes);
  	  // Loop on all the interfaces in my partition
  	  for (std::map<std::vector<int>, simpleInterface>::iterator it = interfacesMap.begin(); it != interfacesMap.end(); ++it) {
  	    std::vector<int> nodes1 = it -> first;
  	    simpleInterface & interface1 = it -> second;
  	    // If it's a periodic BC and we haven't found a partner element yet
  	    if ((interface1.getPhysicalTag()==1)&&(interface1.getElement(1)==NULL)){
  	      // Check if they are paired
  	      bool paired = pairPeriodic(meshNodes, nodes1, nodes);
  	      if(paired){
  		interface1._elements[1] = &el;
  		interface1._closureId[1] = j+nsides;
  	      }
  	    } // end if periodic
  	  } // end loop on interfaces
  	} // end try
  	catch(const std::out_of_range &oor){ continue;} // this interface is not in my partition AND in ghostPeriodic
      } // end catch
    }
  }

  //
  // Check to make sure there are no errors. All interfaces should
  // have neighbors on the right and left.
  //
  for (std::map<std::vector<int>, simpleInterface>::iterator it = interfacesMap.begin(); it != interfacesMap.end(); ++it) {
    simpleInterface & interface = it -> second;
    //printf("interface el1=%i and el2=%i\n",interface.getElement(0)->getId(),interface.getElement(1)->getId());
    if((interface.getElement(0)==NULL)||(interface.getElement(1)==NULL)){
      printf("error in interfaces creation !!!\n");
      exit(1);
    }
  }
    
  //
  // Store the interfaces into the interface vector
  //  
  interfaces.clear();
  interfaces.reserve(interfacesMap.size());
  for (std::map<std::vector<int>, simpleInterface>::iterator it = interfacesMap.begin(); it != interfacesMap.end(); ++it) {
    interfaces.push_back(it->second);
  }
}

void simpleMesh::buildElementMap(int elem_type){
  /* Build a map: element ID -> element index in order of the U
     matrix */
  const std::vector<simpleElement> &elements = getElements(elem_type);
  for(int e = 0; e <elements.size(); e++){
    const simpleElement &el = elements[e];
    //printf("e:%i, id:%i\n", e, el.getId());
    _elementMap[el.getId()] = e;
  }
}
				      
void simpleMesh::buildCommunicators(int elem_type){

  /* Build the necessary maps and matrices to communicate between
     processes.

     Create ghostElementMap. If el2 of an interface belongs to another
     partition, store a unique number. This will be used to access
     back columns of U/A in the limiting procedure (using the
     neighbors vector)

     Create ghostInterfaces matrix of size 3 X # ghost interfaces. It
     stores the local interface index, the partition where el2 lives,
     and a unique tag to reference this interface during MPI
     communication.

     Create ghostElementSend matrix containing the element index to
     send to other partitions, the number of that partition to send
     to, and its global id (as tag for communication).

     Create ghostElementRecv matrix containing the element index to
     store an element from another partition, the source partition of
     that element, and its global id (as tag for communication).

  */

  // Count the number of interfaces which have el1 and el2 on
  // different partitions and fill the ghostElementMap
  int count = 0;
  // counts the number of interfaces two partitions have in common
  int* sharedCount = new int[_numprocs];  for(int k=0; k<_numprocs; k++){sharedCount[k]=0;}
  const std::vector<simpleElement> &elements = getElements(elem_type);
  int N_E = elements.size(); // number of elements in my partition
  for(int i = 0; i < _interfaces.size(); i++){
    const simpleInterface &face = _interfaces[i];  // get the interface
    const simpleElement *el1 = face.getElement(0); // get the element on one side
    const simpleElement *el2 = face.getElement(1); // get the element on the other side
    if( _myid != el2->getPartition()){
      _ghostElementMap[el2->getId()] = N_E + count;
      sharedCount[el2->getPartition()] = sharedCount[el2->getPartition()]++;
      count++;
    }
  }
  _N_ghosts = count;
  
  // initialize the vectors
  _ghostInterfaces  = new int[3*_N_ghosts];
  _ghostElementSend = new int[3*_N_ghosts];
  _ghostElementRecv = new int[3*_N_ghosts];

  // Now fill those vectors
  count = 0;
  for(int t = 0; t < _interfaces.size(); t++){
    const simpleInterface &face = _interfaces[t];  // get the interface
    const simpleElement *el1 = face.getElement(0); // get the element on one side
    const simpleElement *el2 = face.getElement(1); // get the element on the other side
    int partition1 = el1->getPartition();
    int id1 = el1->getId();
    int partition2 = el2->getPartition();
    int id2 = el2->getId();
    if( partition1 != partition2){

      _ghostInterfaces[count*3+0] = t;
      _ghostInterfaces[count*3+1] = partition2; // dest/source partition
      _ghostInterfaces[count*3+2] = uniqueTag(id1,id2,sharedCount[partition2]);

      _ghostElementSend[count*3+0] = _elementMap.at(id1);
      _ghostElementSend[count*3+1] = partition2;
      _ghostElementSend[count*3+2] = id1;

      _ghostElementRecv[count*3+0] = _ghostElementMap.at(id2);
      _ghostElementRecv[count*3+1] = partition2;
      _ghostElementRecv[count*3+2] = id2;      
      
      count++;	
    } // end if condition on partitions
  } // end for loop
  delete[] sharedCount;
}
				      

// Return a matrix which is N_E x D containing the element centroids
fullMatrix<scalar> simpleMesh::getElementCentroids(const int N_E, const int ncorners, const fullMatrix<scalar> XYZNodes){

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

  if      (typeElement=="lin"){ cartesian = true;}
  else if (typeElement=="qua"){
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

  // set the local bool
  _cartesian = cartesian;
  
  return cartesian;
}

void simpleMesh::setDx(const int N_N, const int N_E, const fullMatrix<scalar> &XYZCen, const fullMatrix<scalar> &XYZNodes){
  /* This function finds the minium Dx of all the elements in the
     mesh. It is later used for the CFL condition and the adaptive
     time-stepping.

     For 1D, it's straight forward.  For 2D, it calculates the minimum
     distance from the cell center to edges of the element and it sets
     Dx to twice that value. This means Dx will be equal to the
     conventional Dx for a square element and it will be equal to the
     diameter of the inscribed circle for a triangular element.

  */

  _Dx = 0;
  
#ifdef ONED
  _Dx = XYZNodes(1,0) - XYZNodes(0,0); // initialize first delta x
  // loop on all elements to find the miniumum
  for(int e=0; e<N_E; e++){
    scalar dx = XYZNodes(1,e) - XYZNodes(0,e);
    if (_Dx > dx) _Dx = dx;    
  }

#elif TWOD

  // Distance from first element center to first edge
  scalar d = distance_to_edge(XYZCen(0,0),XYZCen(0,1),XYZNodes(0,0*D+0),XYZNodes(0,0*D+1),XYZNodes(1,0*D+0),XYZNodes(1,0*D+1));
  _Dx = d;

  // Find the other distance in the mesh
  for(int e=0; e<N_E; e++){           // Loop on all the elements
    for(int nn=0; nn<N_N; nn++){     // Loop on all the edges
      d = distance_to_edge(XYZCen(e,0),XYZCen(e,1),XYZNodes(nn,e*D+0),XYZNodes(nn,e*D+1),XYZNodes((nn+1)%N_N,e*D+0),XYZNodes((nn+1)%N_N,e*D+1));
      if(_Dx > d) _Dx = d;
    }
  }

  // Multiply by 2 to get full length of cell
  _Dx = 2*_Dx;
  
#elif THRD
  // distance from point to surface?
#endif
}

inline scalar simpleMesh::distance_to_edge(scalar x0,scalar y0,scalar x1,scalar y1,scalar x2,scalar y2){
  // Distance from center to an edge (from http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html)
  return fabs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}


simpleMesh::~simpleMesh(){
  //delete[] _boundary;
  //delete _boundaryIdx;
}

