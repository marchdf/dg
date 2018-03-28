/*!
  \file simpleMesh.cc
  \brief Function definitions for the simpleMesh class.
  \brief Class deals with the mesh.
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include "simpleMesh.h"
#include "GmshDefines.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <algorithm>
#include "polynomialBasis.h"
#include "Gauss.h"
#include <stdexcept>      // std::out_of_range
#include <bitset>

bool pairPeriodic(const fullMatrix<double> &meshNodes, std::vector<int> nodes1, std::vector<int> nodes2){

  /*!
    \brief Figure out if interfaces are periodic pairs
    \param[in] meshNodes coordinates of all the nodes in the mesh
    \param[in] nodes1 node coordinates of the first interface
    \param[in] nodes2 node coordinates of the second interface
    \return true if periodic pairs
    \section Description
    Given a list of nodes on an interface and their coordinates,
    figure out if they are periodic pairs.

    Used by BuildInterfaces.
  */
  int verbose = 0;
  bool paired = false;
  double eps = 1e-8;
  
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
  if (verbose == 1)
    {
      printf("(xs1,ys1)=(%f,%f) and (xe1,ye1)=(%f,%f)\n",xs1,ys1,xe1,ye1);
      printf("(xs2,ys2)=(%f,%f) and (xe2,ye2)=(%f,%f)\n",xs2,ys2,xe2,ye2);
    }
  // Make sure not to find myself
  if ((fabs(xs1-xs2)<eps)&&(fabs(xe1-xe2)<eps)&&(fabs(ys1-ys2)<eps)&&(fabs(ye1-ye2)<eps)){ return paired;}
  // vertical interface
  else if (fabs(ys1-ye1)>eps){
    // if their start and end nodes match, pair them
    if ((fabs(ys1-ys2)<eps)&&(fabs(ye1-ye2)<eps)){ 
      paired = true;
      if (verbose == 1){printf("Achieved periodic vertical interface\n"); } }
    }
  // horizontal interface
  else if (fabs(xs1-xe1)>eps){
    // if their start and end nodes match, pair them
    if ((fabs(xs1-xs2)<eps)&&(fabs(xe1-xe2)<eps)){ 
      paired = true;
      if (verbose == 1){printf("Achieved periodic horizontal interface\n"); } }
  }
#endif
#ifdef THREED
  //printf("In PariPeriodic, 3D section\n");
  //Copying Marc's strategy from 2D section.
  //For each of the two node lists, find min and max x,y,z coordinates in each list.
  //If these six items match, then the interface is an interior interface.
  //If there is a mismatch, the edge is periodic
  scalar xs1 = meshNodes(0,nodes1[0]);
  scalar xe1 = meshNodes(0,nodes1[0]);
  scalar ys1 = meshNodes(1,nodes1[0]);
  scalar ye1 = meshNodes(1,nodes1[0]);
  scalar zs1 = meshNodes(2,nodes1[0]);
  scalar ze1 = meshNodes(2,nodes1[0]);
  scalar xs2 = meshNodes(0,nodes2[0]);
  scalar xe2 = meshNodes(0,nodes2[0]);
  scalar ys2 = meshNodes(1,nodes2[0]);
  scalar ye2 = meshNodes(1,nodes2[0]);
  scalar zs2 = meshNodes(2,nodes2[0]);
  scalar ze2 = meshNodes(2,nodes2[0]);
  for(int j=1; j<nodes1.size(); j++){
    //printf("j=%d\t",j);
    xs1 = MIN(xs1,meshNodes(0,nodes1[j]));
    xe1 = MAX(xe1,meshNodes(0,nodes1[j]));
    ys1 = MIN(ys1,meshNodes(1,nodes1[j]));
    ye1 = MAX(ye1,meshNodes(1,nodes1[j]));
    zs1 = MIN(zs1,meshNodes(2,nodes1[j]));
    ze1 = MAX(ze1,meshNodes(2,nodes1[j]));
    xs2 = MIN(xs2,meshNodes(0,nodes2[j]));
    xe2 = MAX(xe2,meshNodes(0,nodes2[j]));
    ys2 = MIN(ys2,meshNodes(1,nodes2[j]));
    ye2 = MAX(ye2,meshNodes(1,nodes2[j]));
    zs2 = MIN(zs2,meshNodes(2,nodes2[j]));
    ze2 = MAX(ze2,meshNodes(2,nodes2[j]));
  }
  //09/18/2017: To prevent false positives when pairing periodic faces,
  //I'm going to calculate the dominant face-normal direction. If 
  //the periodicity direction doesn't match the face-normal direction,
  //then it can't be a periodic interface
  int facenorm_dom = -1;
  scalar dx = fabs(xe1-xs1);
  scalar dy = fabs(ye1-ys1);
  scalar dz = fabs(ze1-zs1);
  scalar eps_dom = pow(10,-8);
  //I'm not sure that this arrangement is reliable,
  //but it is working well for now, maybe I'll come
  //back when a test problem breaks it.
  if (dx > eps_dom && dy > eps_dom)
    {
      facenorm_dom = 2; //can be z-periodic, no other options
    }
  else if (dx > eps_dom && dz > eps_dom)
    {
      facenorm_dom = 1; //can be y-periodic, no other options
    }
  else if (dy > eps_dom && dz > eps_dom)
    {
      facenorm_dom = 0; //can be x-periodic, no other options
    }
  else
    {
      printf("In pairPeriodic, the boundary interface does not have a dominant direction\n");
    }
  
  //printf("\n");
  //Maxima and minima have been established. Now, check for coincidence
  if ((fabs(xs1-xs2)<eps)&&(fabs(xe1-xe2)<eps) && (fabs(ys1-ys2)<eps)&&(fabs(ye1-ye2)<eps) && (fabs(zs1-zs2)<eps)&&(fabs(ze1-ze2)<eps))
    {
      if (verbose>0){    printf("Interface locations coincide perfectly, this is an interface against itself\n");}
      return paired;
    }
  else if (fabs(xs1-xs2)>eps && fabs(xe1-xe2)>eps && facenorm_dom == 0)
    {
      //x coordinates mismatch; if y and z coordinates match, it's an x-periodic interface
      if ((fabs(ys1-ys2)<eps) && (fabs(ye1-ye2)<eps) && (fabs(ze1-ze2)<eps) && (fabs(zs1-zs2)<eps))
	{
	  //Here, it is still possible that the faces lie beside eachother, like
	  //a couple of tiles on the kithchen floor
	  if (fabs(xs1-xe2)<eps || fabs(xe1-xs2)<eps)
	    {
	      //This means they share an edge. NOT a periodic pair
	      if (verbose > 0) {printf("The faces share an edge, not x-periodic\n");}
	    }
	  else
	    {
	      paired = true;
	      if (verbose > 0) {printf("\tx-periodic\t");}
	    }
	  /*
	  paired = true;
	  printf("\tx-periodic\t");
	  */
	  }
      
    }
  else if (fabs(ys1-ys2)>eps && fabs(ye1-ye2)>eps && facenorm_dom == 1)
    {
      //y coordinates mismatch; if x and z coordinates match, it's a possible y-periodic interface
      if ((fabs(xs1-xs2)<eps) && (fabs(xe1-xe2)<eps) && (fabs(ze1-ze2)<eps) && (fabs(zs1-zs2)<eps))
	{
	  //Here, it is still possible that the faces lie beside eachother, like
	  //a couple of tiles on the kithchen floor
	  if (fabs(ys1-ye2)<eps || fabs(ye1-ys2)<eps)
	    {
	      //This means they share an edge. NOT a periodic pair
	      if (verbose > 0) {printf("The faces share an edge, not y-periodic\n");}
	    }
	  else
	    {
	      paired = true;
	      if (verbose > 0) {printf("\ty-periodic\t");}
	    }
	}

    }
  else if (fabs(zs1-zs2)>eps && fabs(ze1-ze2)>eps && facenorm_dom == 2)
    {
      //z coordinates mismatch; if x and y coordinates match, it's a z-periodic interface
      if ((fabs(xs1-xs2)<eps) && (fabs(xe1-xe2)<eps) && (fabs(ye1-ye2)<eps) && (fabs(ys1-ys2)<eps))
	{
	  //Here, it is still possible that the faces lie beside eachother, like
	  //a couple of tiles on the kithchen floor
	  if (fabs(zs1-ze2)<eps || fabs(ze1-zs2)<eps)
	    {
	      //This means they share an edge. NOT a periodic pair
	      if (verbose > 0) {printf("The faces share an edge, not z-periodic\n");}
	    }
	  else
	    {
	      paired = true;
	      if (verbose > 0) {printf("\tz-periodic\t");}
	    }
	  /*
	  paired = true;
	  printf("\tz-periodic\t");
	  */
	}

    }
  /*
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
  //front/back interface
  else if (fabs(zs1-ze1)>eps){
    // if their start and end nodes match, pair them
    if ((fabs(zs1-zs2)<eps)&&(fabs(ze1-ze2)<eps)){ paired = true;}
  }
  */
#endif

  return paired;  
}


int uniqueTag(int id1, int id2, int N=0){
  /*!
    \brief Generates a unique tag given two numbers
    \param[in] id1 first number
    \param[in] id2 second number
    \param[in] N total number of pairs

    \section Description
    This is not used anymore... but I keep it because it's kind of cool

    Given two numbers (the global id of el1 and el2 of an interface),
    create a unique number. Ideally I would use a perfect hash. But
    for that I would need a 64-bit integer to represent all pairs of
    32-bit integers. For a good discussion on this see
    http://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way
    http://stackoverflow.com/questions/682438/hash-function-providing-unique-uint-from-an-integer-coordinate-pair
    http://stackoverflow.com/questions/11786635/fast-bi-directional-hash-of-two-integers-in-c
  */
  
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
  /*!
    \brief Cantor pairing function.
    \param[in] k1 first number
    \param[in] k2 second number
    \section Description
    Generates a number based on two other ones.
    
    This will overflow when for (k1,k2)=(65535, 65535). You can
    represent all 16-bit integers using 32-bit integers but you can't
    go further than that without using a 64-bit integer.
  */
  return 0.5*(k1+k2)*(k1+k2+1)+k2;
}
    
void simpleMesh::load (const char *fileName)
{
  /*!
    \brief Load a mesh file
    \param[in] fileName mesh file to load
  */
  int verbose = 0;
  std::ifstream input (fileName);
  std::string line;
  getline (input, line);
  getline (input, line);
  getline (input, line);
  getline (input, line);
  if (line!="$Nodes")
    printf("invalid file format, did not find Nodes line\n");
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
    printf("invalid file format, did not find EndNodes line\n");
  getline (input, line);
  if (line!="$Elements")
    printf("invalid file format, did not find Elements line\n");
  int nelements;
  input >> nelements;
  std::vector<int> enodes;
  _elements.resize(MSH_NUM_TYPE);      // holds elements in my partition
  _otherElements.resize(MSH_NUM_TYPE); // holds elements in other partitions
  getline (input, line);
  if (_myid == 0) printf("In load subroutine: nelements=%d\n", nelements);
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
      if(_myid == 0 && verbose > 0){printf("Id=%i: j=%i, tag=%i, ptag=%i, partition=%i \n",elementId, j, tag,ptag, partition);}
    }
    int idNode;
    while (sline >> idNode) {
      enodes.push_back(idNode-1);
    }
    if (verbose >0)
      {
	printf("The local enodes structure:\n");
	for (int i = 0; i < enodes.size(); i++)
	  {
	    printf("\tenodes[%d]=%d\n", i, enodes[i]);
	  }
      }
    // Exit if the partition number is larger than the total number of processors
    if(_numprocs < partition)
      {
	printf("Not enough processors for the mesh partitions. Exiting\n"); 
	exit(1);
      }
    // Only store the element if it's in my partition
    if(_myid == partition-1) 
      {
	_elements[elementType].push_back (simpleElement(elementId, ptag, partition-1, enodes));
	if (verbose > 0) {printf("Pushed a new element: elementType=%d, elementId=%d, physical=%d, partition=%d. elements size now = %ld\n", elementType, elementId, ptag, partition-1, _elements.size());}
      }
    // Otherwise store the element in the otherElements
    else 
      {
	_otherElements[elementType].push_back (simpleElement(elementId, ptag, partition-1, enodes));
	if (verbose > 0) {printf("Pushed a new OtherElement: elementType=%d, elementId=%d, physical=%d, partition=%d, otherElements size now = %ld\n", elementType, elementId, ptag, partition-1, _otherElements.size());}
      }
  }
  getline (input, line);
  if (line!="$EndElements")
    printf("invalid file format, could not find EndElements line\n");
  /*
  if (_myid == 0)
    {
      printf("Inside load routine, myid == 0: Here is the _elements structure, type 3:\n");
      for (int e = 0; e < _elements.size(); e++)
	{
	  int elementType = 3;
	  printf("elementId=%d, physical=%d, partition=%d\n", _elements[elementType][e].getId(), _elements[elementType][e].getPhysicalTag(), _elements[elementType][e].getPartition());
	}
    }
  
  if (_myid == 1)
    {
      printf("Inside load routine, myid == 1: Here is the _elements structure, type 3:\n");
      for (int e = 0; e < _elements.size(); e++)
	{
	  int elementType = 3;
	  printf("elementId=%d, physical=%d, partition=%d\n", _elements[elementType][e].getId(), _elements[elementType][e].getPhysicalTag(), _elements[elementType][e].getPartition());
	}
    }
  */
}

void simpleMesh::buildNormals (int typeInterface, int typeElement, int PeriCo2D, std::vector<std::vector<int> > elemNodes, int order, int N_N)
{
  /*!
    \brief Build the normals to all the elements
    \param[in] typeInterface the associated key number referencing the type of interface
    \param[in] typeElement the associated key number referencing the type of element
    \param[in] PeriCo2D whether or not to execute periodicity fix on node locations
    \section Description
    Build the normals to all the elements
  */
  
  //
  // Initialize the elements, interfaces, nodal functions and closures
  //
  int verbose = 0;
  if (verbose>0 && _myid == 0){
    printf("Entered buildNormals subroutine\n");
    printf("The _nodes structure:\n");
    for (int j = 0; j < _nodes.size2(); j++)
      {
	for (int i = 0; i < _nodes.size1(); i++)
	  {
	    printf("_nodes(a=%d, index=%d) = %f, ",i,j,_nodes(i,j));
	  }
	printf("\n");
      }
  }
  
  const std::vector<simpleElement> &elements = getElements(typeElement);
  const std::vector<simpleInterface> &interfaces = getInterfaces();
  const polynomialBasis &basis  = *polynomialBases::find (typeElement);  // for the element
  const std::vector<std::vector<int> > &closures = basis.closures;

  // Define some numbers for clarity
  if (verbose > 0) {printf("Attempting to fetch N_s from first element's information:\n"); fflush(stdout);}
  int N_s = elements[0].getNbNodes(); // number of nodes on an element          (i index)
  if (verbose > 0) {printf("Fetched N_s = %d\n", N_s); fflush(stdout);} 
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
  for(int g = 0; g < N_s; g++){
    basis.df(points(g,0),points(g,1),points(g,2),grads);
    for(int alpha = 0; alpha < D; alpha ++){
      for(int i = 0; i < N_s; i++){
   	dphi(g*D+alpha,i) = grads[i][alpha];  // see paper for indexing p.6
      }	  
    }    
  }
  if (verbose > 0)
    {

  printf("In buildNormals: Basis in the reference element:\n");
  for (int j = 0; j < N_s; j++)
    {
      printf("node %d: (xi,eta,zeta) = (%f, %f, %f), nonzero basis functions are\n", j, points(j,0), points(j,1), points(j,2));
      for (int k = 0; k < N_s; k++)
	{
	  if (fabs(phi(j,k) > pow(10,-10)))
	    {
	      printf("\tphi[index=%d] = %f\n",k, phi(j,k));
	    }
	} 
    }
    }
  //
  // Calculate the normals
  //
  _normals.resize(D,N_I);
  _normals.setAll(0.0);
  double* n = new double[D];
  
   //PEJ 11/06/2017: Instead of working with _nodes,
    //this code will work with MeshNodes, which can perform periodicity
    //fix on _nodes
  //ALSO NEED TO MAKE THE SECOND PERIFIX ADJUTMENT, REQUIRES ELEMENTS TO BE LOADED (done above)
  fullMatrix<double> MeshNodes = PeriFixNodes(PeriCo2D, _nodes);
  MeshNodes = PeriFixNodes_HiOMesh(PeriCo2D, MeshNodes, elemNodes, typeElement, order, N_N);
  if (verbose>0 && _myid == 0){
    printf("In buildnormals: The _nodes structure:\n");
    for (int j = 0; j < _nodes.size2(); j++)
      {
	for (int i = 0; i < _nodes.size1(); i++)
	  {
	    printf("_nodes(a=%d, index=%d) = %f, ",i,j,_nodes(i,j));
	    }
	  printf("\n");
	}
      printf("In buildnormals: The MeshNodes structure after both periodicity corrections:\n");
      for (int j = 0; j < MeshNodes.size2(); j++)
	{
	  for (int i = 0; i < MeshNodes.size1(); i++)
	    {
	      printf("MeshNodes(a=%d, index=%d) = %f, ",i,j,MeshNodes(i,j));
	    }
	  printf("\n");
	}
    }

  // Loop on all the interfaces
  for(int i = 0; i < N_I; i++){
    if (verbose == 1)
      {
	printf("Entered interface loop to determing normal, i=%d out of %d\n",i,N_I);
      }
    // Get the interface
    const simpleInterface &face = interfaces[i];

    // Get some information: the element, the closure.
    //first element of the interface is always
    //on the partition, so no parallelization worries here.
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
    	//XYZNodes(j, alpha) = _nodes(alpha, inode);
	XYZNodes(j, alpha) = MeshNodes(alpha, inode);
      }
    }

    if (verbose > 0 && D == 3)
      {
	printf("The interface's el0 XYZNodes:\n");
	for (int j = 0; j < N_s; j++)
	  {
	    printf("node %d: (x,y,z) = (%f, %f, %f)\n", j, XYZNodes(j,0), XYZNodes(j,1), XYZNodes(j,2));
	  }
      }

    // Jacobian matrix of the element
    fullMatrix<double> Jac(D,N_s*D);
    Jac.gemm(XYZNodes.transpose(),dphi.transpose());
    for (int j = 0; j < N_s; j++) //node where we are calculating jacobian
      {
	for (int row = 0; row < D; row++) //row of the jacobian
	  {
	    for (int col = 0; col < D; col++) //column of the jacobian
	      {
		//we seek d(x_row)/d(xi_col)
		scalar sum = 0;
		for (int index = 0; index < N_s; index++) //influential basis index
		  {
		    sum += dphi(j*D+col , index) * XYZNodes(index, row);
		  }
		Jac(row, j*D + col) = sum;
	      }
	  }
      }
    //PEJ: I think this is the jacobian at each solution node,
    //meaning it accounts for variation across the element.
    //Could be useful in future

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
#ifdef THREED
      //For inverse Jacobian, try straight fullMatrix inversion
      //Take Jacobian, perform inverse routine, store result in invJac:
      fullMatrix<double> JacNode(D,D);
      fullMatrix<double> invJacNode(D,D);
      for (int a1 = 0; a1 < D; a1++)
	{
	  for (int a2 = 0; a2 < D; a2++)
	    {
	      //Not sure if this indexing is right
	      JacNode(a1,a2) = Jac(a1, k*D + a2);
	    }
	}
      //Invert JacNode, store result in invJacNode
      JacNode.invert(invJacNode);
      //Now, relay the result to global storage. Marc Jumbles the organization
      //in 2D case, so I will repeat for 3D
      for (int a1 = 0; a1 < D; a1++)
	{
	  for (int a2 = 0; a2 < D; a2++)
	    {
	      //not sure if this indexing is right
	      invJac(a2, k*D+a1) = invJacNode(a1,a2);
	    }
	}
      if (verbose == 1)
	{
	  printf("In BuildNormals: the Jacobian matrix at node %d:\n",k);
	  for (int row = 0; row < D; row++)
	    {
	      printf("row %d = %f\t%f\t%f\n", row, JacNode(row,0), JacNode(row,1), JacNode(row,2));
	    }
	  printf("In BuildNormals: the inverse Jacobian matrix at node %d:\n",k);
	  for (int row = 0; row < D; row++)
	    {
	      printf("row %d = %f\t%f\t%f\n", row, invJacNode(row,0), invJacNode(row,1), invJacNode(row,2));
	    }
	}
    
#endif
    } //end k loop

    // Loop on the nodes of the interface and calculate the normal
    for(int k = 0; k < cl.size(); k++){
      for(int alpha = 0; alpha < D; alpha++){
	for(int a = 0; a < D; a++){ // loop for matrix-matrix mult of invJac and dphi
	  n[alpha] += invJac(a,cl[k]*D+alpha)*dphi(cl[k]*D+a,cl[k]);
	}
      }
    }

    if (verbose == 1)
      {
	printf("Calculated interface normal, pre-normalization:\t");
	printf("normal = (");
	for (int a = 0; a < D; a++)
	  {
	    printf("%f,  ",n[a]);
	  }
	printf(")\n");
      }

    // Normalize and store
    scalar norm2 = 0.0;
    for(int alpha=0; alpha < D; alpha++){ norm2 += n[alpha]*n[alpha];}
    for(int alpha=0; alpha < D; alpha++){
      _normals(alpha,i) = (scalar)n[alpha]/sqrt(norm2);
    }
    if (verbose == 1)
      {
	printf("Calculated interface normal, post-normalization:\t");
	printf("normal = (");
	for (int a = 0; a < D; a++)
	  {
	    printf("%f,  ",_normals(a,i));
	  }
	printf(")\n");
      }
  } // end loop on interfaces
  delete[] n;
}


void simpleMesh::writeSolution (const scalar* solution, const int N_s, const int N_E, int type, std::vector<std::string> fnames, std::vector<std::string> names, int step, double time) const
{
  /*!
    \brief Write the solution to a file
    \param[in] solution array of solution to output
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] type element type to output
    \param[in] fnames vector of output file names (w/o extension)
    \param[in] names vector of output field names
    \param[in] step time step number
    \param[in] time time value
    \section Description
    Pretty obvious what this does.
  */

  const std::vector<simpleElement> &list = _elements[type];
  if (list.size() != N_E)
    printf("bad solution for this element\n");

  // Loop on number of fields for output
  int _N_F = fnames.size();
  for(int fc=0; fc< _N_F; fc++){
    std::ofstream output;
    std::string filename = fnames[fc];
    char stepstr[21];
    sprintf(stepstr, "%010i", step);
    filename += stepstr;
    filename += ".pos"; // add extension
#ifdef USE_MPI
    char myidstr[21]; // enough to hold all numbers up to 64-bits
    sprintf(myidstr, "%d", _myid);
    filename += myidstr;
#endif
    output.open(filename.c_str()); //output.open(filename.c_str(),std::ios_base::app); 
    output.precision(20);
    output << "$MeshFormat\n2.1 0 8\n$EndMeshFormat\n";
    output << "$ElementNodeData\n";
    output << "1\n\"" << names[fc] << "\"\n";
    output << "1\n" << time << "\n";
    output << "4\n" << step<< "\n1\n" << list.size() << "\n" << _myid << "\n";
    for (int e = 0; e < list.size(); e++) {
      const simpleElement &element = list[e];
      output << element.getId() << " " << element.getNbNodes() <<" ";
      if (element.getNbNodes() != N_s)
      	printf("bad solution for this element\n");
      for (int i = 0; i < element.getNbNodes(); i++) {
	output << solution[(e*_N_F+fc)*N_s+i]<< " ";
      }
      output << "\n";
    }
    output << "$EndElementNodeData\n";
    output.close();
  } // end loop on fields
}

void simpleMesh::readSolution (const int N_s, const int N_E, int type, std::vector<std::string> fnames, std::vector<std::string> names, const int step, double &time, scalar* solution){
  /*!
    \brief Read the solution from input files
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] type element type to output
    \param[in] fnames vector of output file names (w/o extension)
    \param[in] names vector of output field names
    \param[in] step time step number (of file to read)
    \param[out] time time value
    \param[out] solution array of solution to output
    \section Description
    Pretty obvious what this does.
  */

  const std::vector<simpleElement> &list = _elements[type];
  if (list.size() != N_E)
    printf("bad solution for this element\n");

  // Loop on number of fields for output
  int _N_F = fnames.size();
  for(int fc=0; fc< _N_F; fc++){
    std::ifstream input;
    std::string filename = fnames[fc];
    char stepstr[21];
    sprintf(stepstr, "%010i", step);
    filename += stepstr;
    filename += ".pos"; // add extension
#ifdef USE_MPI
    char myidstr[21]; // enough to hold all numbers up to 64-bits
    sprintf(myidstr, "%d", _myid);
    filename += myidstr;
#endif
    input.open(filename.c_str(),std::ifstream::in);
    if(input.is_open()==0){
      std::cout << "No file named " << filename << ". Exiting." << std::endl;
      exit(1);
    }
    std::string line;
    getline(input,line); // ignore $MeshFormat
    getline(input,line); // ignore 2.1 0 8
    getline(input,line); // ignore $EndMeshFormat
    getline(input,line); // ignore $ElementNodeData
    getline(input,line); // ignore 1
    getline(input,line); // ignore name
    getline(input,line); // ignore 1
    input >> time; getline(input,line); // get time
    getline(input,line); // ignore 4
    getline(input,line); // ignore step
    getline(input,line); // ignore 1
    int fsize; input >> fsize; getline(input,line); // test list size
    if(fsize != list.size()){
      std::cout << "File has a different number of elements than current mesh. Exiting." << std::endl;
      exit(1);
    }
    int fmyid; input >> fmyid; getline(input,line); // test partition
    if(fmyid != _myid){
      std::cout << "File has a different partition number than me. Exiting." << std::endl;
      exit(1);
    }

    // Read the data
    for (int e = 0; e < list.size(); e++) {
      const simpleElement &element = list[e];
      int f_elementId;
      int f_NbNodes;
      input >> f_elementId;
      input >> f_NbNodes;
      if ((f_elementId != element.getId()) || (f_NbNodes != element.getNbNodes())){
	std::cout << "File has a problem (element id or element number of nodes). Exiting." << std::endl;
	exit(1);
      }

      // loop on nodal data
      for (int i = 0; i < element.getNbNodes(); i++) {
	input >> solution[(e*_N_F+fc)*N_s+i];
      }
      getline(input,line); // get the rest of the line
    }
    input.close();
  } // end loop on fields
}

void simpleMesh::readSolutionNoTime (const int N_s, const int N_E, int type, std::vector<std::string> fnames, std::vector<std::string> names, const int step, scalar* solution){
  /*!
    \brief Read the solution from input files
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] type element type to output
    \param[in] fnames vector of output file names (w/o extension)
    \param[in] names vector of output field names
    \param[in] step time step number (of file to read)
    \param[out] solution array of solution to output
    \section Description
    Pretty obvious what this does.
  */

  const std::vector<simpleElement> &list = _elements[type];
  if (list.size() != N_E)
    printf("bad solution for this element\n");

  // Loop on number of fields for output
  int _N_F = fnames.size();
  for(int fc=0; fc< _N_F; fc++){
    std::ifstream input;
    std::string filename = fnames[fc];
    char stepstr[21];
    sprintf(stepstr, "%010i", step);
    filename += stepstr;
    filename += ".pos"; // add extension
#ifdef USE_MPI
    char myidstr[21]; // enough to hold all numbers up to 64-bits
    sprintf(myidstr, "%d", _myid);
    filename += myidstr;
#endif
    input.open(filename.c_str(),std::ifstream::in);
    if(input.is_open()==0){
      std::cout << "No file named " << filename << ". Exiting." << std::endl;
      exit(1);
    }
    std::string line;
    getline(input,line); // ignore $MeshFormat
    getline(input,line); // ignore 2.1 0 8
    getline(input,line); // ignore $EndMeshFormat
    getline(input,line); // ignore $ElementNodeData
    getline(input,line); // ignore 1
    getline(input,line); // ignore name
    getline(input,line); // ignore 1
    scalar time; //don't care about time
    input >> time; getline(input,line); // get time
    getline(input,line); // ignore 4
    getline(input,line); // ignore step
    getline(input,line); // ignore 1
    int fsize; input >> fsize; getline(input,line); // test list size
    if(fsize != list.size()){
      std::cout << "File has a different number of elements than current mesh. Exiting." << std::endl;
      exit(1);
    }
    int fmyid; input >> fmyid; getline(input,line); // test partition
    if(fmyid != _myid){
      std::cout << "File has a different partition number than me. Exiting." << std::endl;
      exit(1);
    }

    // Read the data
    for (int e = 0; e < list.size(); e++) {
      const simpleElement &element = list[e];
      int f_elementId;
      int f_NbNodes;
      input >> f_elementId;
      input >> f_NbNodes;
      if ((f_elementId != element.getId()) || (f_NbNodes != element.getNbNodes())){
	std::cout << "File has a problem (element id or element number of nodes). Exiting." << std::endl;
	exit(1);
      }

      // loop on nodal data
      for (int i = 0; i < element.getNbNodes(); i++) {
	input >> solution[(e*_N_F+fc)*N_s+i];
      }
      getline(input,line); // get the rest of the line
    }
    input.close();
  } // end loop on fields
}

void simpleMesh::buildBoundary(){
  /*!
    \brief Build a list of special boundaries
    
    \section Description
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
  int verbose = 0;
  //
  // Sort each interface according to type
  //
  std::vector<simpleInterface> rflctive;
  std::vector<simpleInterface> otherone;
  std::vector<simpleInterface> noslip; //PEj insert, 01/20/2016
  std::vector<simpleInterface> nograd; //PEJ inset, 01/20/2016: This is Marc's farfield boundary interfaces
  std::vector<simpleInterface> Anflw; //PEJ 01/21/2016: an inflow condition
  std::vector<simpleInterface> KJet; //PEJ 10/10/2017: The Kushner jet (equal-pressure, low density)
  std::vector<simpleInterface> SubOut; //PEJ 10/10/2017: subsonic outflow
  std::vector<simpleInterface> Homo; //PEJ 10/10/2017: Homogeneous (U=2) boundary condition

  std::vector<int> rflctiveFaceNumber; // these hold the index to the face in the general matrix
  std::vector<int> otheroneFaceNumber; 
  std::vector<int> noslipFaceNumber; //PEJ insert, 01/20/2016
  std::vector<int> nogradFaceNumber; //PEJ, 01/20/2016
  std::vector<int> AnflwFaceNumber; //PEJ, 01/21/2016
  std::vector<int> KJetFaceNumber; //PEJ, 10/01/2017
  std::vector<int> SubOutFaceNumber; //PEJ, 10/01/2017
  std::vector<int> HomoFaceNumber; //PEJ, 10/01/2017

  if (verbose > 0) {printf("Looping over ALL interfaces in buildBoundary subroutine:\n");}
  for(int i = 0; i < N_I; i++)
    {
      //printf("Working interface %d:    ",i);
      const simpleInterface &face = _interfaces[i];
      int physical = face.getPhysicalTag();
      //printf("physical=%d,   ",physical);
      if (3==physical)
	{ 
	  rflctive.push_back(face); 
	  rflctiveFaceNumber.push_back(i);
	  //printf("stacked on reflective vector");
	}
      else if(5==physical)
	{ 
	  otherone.push_back(face); 
	  otheroneFaceNumber.push_back(i);
	  //printf("stacked on otherone vector");
	}
      //PEJ Edit, 01/20/2016:
      else if(7==physical) 
	{ 
	  noslip.push_back(face);   
	  noslipFaceNumber.push_back(i);
	  //printf("stacked on noslip vector");
	}
      else if(9==physical) 
	{ 
	  Anflw.push_back(face);   
	  AnflwFaceNumber.push_back(i);
	  //printf("stacked on Anflw vector");
	}
      else if (2==physical) 
	{
	  nograd.push_back(face);   
	  nogradFaceNumber.push_back(i);
	  //printf("stacked on nograd vector");
	}
      else if (11==physical)
	{
	  KJet.push_back(face);
	  KJetFaceNumber.push_back(i);
	  //	  printf("stacked on KJet vector");
	}
      else if (13==physical)
	{
	  SubOut.push_back(face);
	  SubOutFaceNumber.push_back(i);
	  //	  printf("stacked on KJet vector");
	}
      else if (101==physical)
	{
	  Homo.push_back(face);
	  HomoFaceNumber.push_back(i);
	  //	  printf("stacked on Homogeneous BC vector\n");
	}
      //printf("\n");
    }

  //
  // Number of boundary faces
  //
  //_N_B = rflctive.size() + otherone.size();
  //_N_B = rflctive.size() + otherone.size() + noslip.size();
  //_N_B = rflctive.size() + otherone.size() + noslip.size() + nograd.size();
  //_N_B = rflctive.size() + otherone.size() + noslip.size() + nograd.size() + Anflw.size();
  //_N_B = rflctive.size() + otherone.size() + noslip.size() + nograd.size() + Anflw.size() + KJet.size();
  //_N_B = rflctive.size() + otherone.size() + noslip.size() + nograd.size() + Anflw.size() + KJet.size() + Homo.size();
  _N_B = rflctive.size() + otherone.size() + noslip.size() + nograd.size() + Anflw.size() + KJet.size() + Homo.size() + SubOut.size();
  /*
  printf("reflective boundaries:%d\n",int(rflctive.size()));
  printf("otherone boundaries:%d\n",int(otherone.size()));
  printf("noslip boundaries:%d\n",int(noslip.size()));
  */

  //
  // Match the boundaries
  // 
  _boundary = new int[_N_B];
  //PEJ note: I'm thinking that _boundaryIdx size should match the number of 
  //boundary treatments being used. For example, if some interfaces are no-slip, some are reflective,
  //and some are zero-gradient, the declaration might need to be new int[3]. 
  //_boundaryIdx = new int[1];
  //I'm thinking that _boundaryIdx should be vector<int>, so we can resize it automatically based
  //on how many types of BC are present on the mesh. That's a problem for future Phil

  _boundaryIdx = new int[7]; //7 entries, for 7 types of BC: Change this number if you want more types of BC
  int t = 0; //I think this is a master index for the boundary interface, used to map to actual interface index

  // Reflective BC
  for(int t1=0; t1<rflctive.size(); t1++){
    //printf("Entered reflective BC loop, t=%d\n",t);
    // get the face
    const simpleInterface &face1 = rflctive[t1];
    int idx1 = rflctiveFaceNumber[t1];
    _boundary[t] = idx1;
    t++;
  }
  _boundaryIdx[0] = t; //this is the end point of the reflective interfaces
  
  // Other BC (not used for simulations, just a model to show newcomers how to deal with boundaries)
  for(int t1=0; t1<otherone.size(); t1++){
    //printf("Entered otherone BC loop, t=%d\n",t);
    // get the face
    const simpleInterface &face1 = otherone[t1];
    int idx1 = otheroneFaceNumber[t1];
    _boundary[t] = idx1;
    t++;
  }
  _boundaryIdx[0] = t;
  
  //PEJ Edit: No-slip interfaces:
  for (int t1 = 0; t1<noslip.size();t1++)
    {
      //printf("Entered noslip BC loop, t=%d\n",t);
      //get the face
      const simpleInterface &face1 = noslip[t1];
      int idx1 = noslipFaceNumber[t1];
      _boundary[t] = idx1;
      t++;
    }
  _boundaryIdx[1] = t; //this is the end point of the no-slip interfaces

  //Now, the nograd interfaces
  for (int t1 = 0; t1<nograd.size();t1++)
    {
      //printf("Entered nograd BC loop, t=%d\n",t);
      //get the face
      const simpleInterface &face1 = nograd[t1];
      int idx1 = nogradFaceNumber[t1];
      _boundary[t] = idx1;
      t++;
    }
  _boundaryIdx[2] = t; //this is the end point of the no-grad interfaces

  //Now, the Anflw interfaces
  for (int t1 = 0; t1<Anflw.size();t1++)
    {
      //printf("Entered Anflw BC loop, t=%d\n",t);
      //get the face
      const simpleInterface &face1 = Anflw[t1];
      int idx1 = AnflwFaceNumber[t1];
      _boundary[t] = idx1;
      t++;
    }
  _boundaryIdx[3] = t; //this is the end point of the A inflow interfaces

  //Now, the KJet interfaces
  for (int t1 = 0; t1<KJet.size();t1++)
    {
      if (verbose > 0) {printf("Entered KJet BC loop, t=%d\n",t);}
      //get the face
      const simpleInterface &face1 = KJet[t1];
      int idx1 = KJetFaceNumber[t1];
      _boundary[t] = idx1;
      t++;
    }
  _boundaryIdx[4] = t; //this is the end point of the KJet interfaces

  //Now, the Homogeneous BC interfaces
  for (int t1 = 0; t1<Homo.size();t1++)
    {
      if (verbose > 0) {printf("Entered Homogeneous BC loop, t=%d\n",t);}
      //get the face
      const simpleInterface &face1 = Homo[t1];
      int idx1 = HomoFaceNumber[t1];
      _boundary[t] = idx1;
      t++;
    }
  _boundaryIdx[5] = t; //this is the end point of the Homogeneous BC interfaces

  //Now, theSubSoninc Outflow BC interfaces
  for (int t1 = 0; t1<SubOut.size();t1++)
    {
      if (verbose > 0) {printf("Entered SubOut BC loop, t=%d\n",t);}
      //get the face
      const simpleInterface &face1 = SubOut[t1];
      int idx1 = SubOutFaceNumber[t1];
      _boundary[t] = idx1;
      t++;
    }
  _boundaryIdx[6] = t; //this is the end point of the Subsonic Outflow BC interfaces

  //END PEJ Edit
  if (verbose > 0) {printf("boundaryIdx[0] = %d, boundaryIdx[1] = %d, boundaryIdx[2] = %d, boundaryIdx[3] = %d, boundaryIdx[4] = %d, boundaryIdx[5] = %d, boundaryIdx[6] = %d\n", _boundaryIdx[0], _boundaryIdx[1], _boundaryIdx[2], _boundaryIdx[3], _boundaryIdx[4], _boundaryIdx[5], _boundaryIdx[6]);}

  if (verbose == 1)
    {
      printf("Finished with build boundary routine: here is the boundary map:\n");
      for (int t = 0; t < _N_B; t++)
	{
	  printf("\tboundary interface %d: global interface address is %d\n",t,_boundary[t]);
	}
    }
 
}


void simpleMesh::buildBoundaryElementShift1D(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes){ 
  /*!
    \brief Build 1D shifts to bring an element next to his neighbor
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] XYZNodes element node coordinates
    \section Description
    Objective: create BoundaryElemShift
    elem1         | elem2         | xshift
    (tgt element)   his neighbor    shifts to bring elem2 next to elem1
  */

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
  /*!
    \brief Build 2D shifts to bring an element next to his neighbor
    \param[in] order DG order
    \param[in] XYZNodesF nodal coordinates of interfaces
    \param[in] ElementMap map from element id to element index
    \section Description
    Objective: create BoundaryElemShift
    elem1         | elem2         | xshift | yshift
    (tgt element)   his neighbor    shifts to bring elem2 next to elem1
  */
  
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

int simpleMesh::SeekGhost(/*const simpleElement* elHome, int Actual_om_on, */int GmshId_on, int GmshId_off)
{
  /*!
    \brief identify the [ghost index of the ghost element] correspinding to the [flesh element with gmsh address GmshId_off]
    \param[in] index in (0..N_E) of the element on the myid partition
    \param[in] GmshId_on the integer value that the gmsh mesh generator assigns to myid partition's neighbor of GmshId_off
    \param[in] GmshId_off the integer value that the gmsh mesh generator assigns to the element we are working with
    \param[out] ghostId the integer index (between N_E and N_E+N_ghosts) that tells myid processor where to store ghost element information
   */
  int ghostId = -1;
  for (int j = 0; j < _ghostElement_Key.size(); j++)
    {
      int key = _ghostElement_Key[j];
      if (key == GmshId_off)
	{
	  //We have found a key value corresponding to the flesh element with gmsh address GmshId_off.
	  //Problem: there may be multiple entries in _ghostElement_Key that match GmshId_off (partition corner intersections).
	  //Secondary check to see that elId=GmshId_on and elIndex=_ghostElement_HNe[j] are same element.
	  //HNe stands for home element, it is the gmsh index of the on-partition neighbor of flesh element w/ gmsh id "key"
	
	  int om_off = _ghostElement_Val[j]; //myid index of the ghost replica of GmshId_off, corresponds to _ghostElement_Key[j]
	  int GmshId_HNe = _ghostElement_HNe[j];
	  //  printf("\nIn SeekGhost: pr=%d, GmshId_on=%d, GmshId_off=%d, and _ghostElement_Key[%d]=%d. Possible Match\n", _myid, GmshId_on, GmshId_off, j, key);
	  if (GmshId_on == GmshId_HNe)
	    {
	      //  printf("\tWe have a match with GmshId_HNe=%d, GmshId_on=%d, om_off=%d\n",GmshId_HNe, GmshId_on, om_off);
	      ghostId = om_off;
	    }
	  else
	    {
	      //   printf("\tYou are NOT THE FATHER! GmshId_HNe=%d, GmshId_on=%d, om_off=%d\n",GmshId_HNe, GmshId_on, om_off);
	      //NO MATCH! try again :)
	    }
	}
    }
  return ghostId;
}

void simpleMesh::buildNeighbors(int N_N, int N_E)
{
  /*!
    \brief Build the neighbors
    \section Description
    Objective: find the neighbors of each element
    get _neighbors, a N_N x N_E vector:
    | neighbor1 | ...
    | neighbor2 | ...
  */

  // Allocate neighbors, set to zero
  int verbose = 0;
  if (verbose == 1) {printf("Entered buildNeighbors: N_E=%d, N_N=%d\n", N_E, N_N);}

  _neighbors = new int[N_N*N_E]; for(int k=0; k < N_N*N_E; k++){_neighbors[k]=0;}

  int N_I = _interfaces.size();       // number of interfaces                   (i index)
  if (verbose == 1) {printf("Entered buildNeighbors: N_E=%d, N_N=%d, N_I=%d\n", N_E, N_N, N_I);}
  double eps = 1e-7; // used for comparisons of smallness
  
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
    else                          {el2num = SeekGhost(/*el1, el1num, */el1->getId(), el2->getId());}
    //else                          {el2num = _ghostElementMap.at(el2->getId());}
    
    if (verbose == 1) 
      {
	printf("Working on interface %d: el1=%d, el2=%d\n", i,el1num,el2num);
      }
    
    // If el1num = el2num, we are at some kind of boundary
    if (el2num == el1num){
      int physical = face.getPhysicalTag();
      // give elnum2 the negative value of the physical. The idea is
      // that when using the neighbors in the limiting procedure, we can
      // know which boundary we are at and use the appropriate nodal
      // values for the reconstruction
      
      //if ((physical == 2)||(physical == 3)||(physical == 4)){el2num = -physical;}
      //PEJ 10/31/2017: Adding some more boundaries that should trigger sensor
      if ((physical == 2)||(physical == 3)||(physical == 4)||(physical == 101)||(physical == 9)){el2num = -physical;}
    }
    
    // Not a cartesian mesh
    if(!_cartesian){
      //  if (verbose == 1) {printf("\n\nMesh is NOT Cartesian\n\n");}
      _neighbors[el1num*N_N+nn[el1num]] = el2num; nn[el1num]++;
      // do the reverse if el2 belongs to me
      if (el2num>=0){ if(el2->getPartition()==_myid){ _neighbors[el2num*N_N+nn[el2num]] = el1num; nn[el2num]++;}}
    }
    else if (_cartesian){ // sort in LRDU order
      //      if (verbose == 1) {printf("\n\nMesh is Cartesian\n\n");}
      double nx = _normals(0,i);
#ifdef ONED
      {
	double ny = 0;
	//printf("e1=%i, e2=%i, nx=%e, ny=%e\n",el1->getId(),el2->getId(),nx,ny);
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
      }
#elif TWOD
      {
	double ny = _normals(1,i);
	//printf("e1=%i, e2=%i, nx=%e, ny=%e\n",el1->getId(),el2->getId(),nx,ny);
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
      }
#elif THREED
      {
	double ny = _normals(1,i);
	double nz = _normals(2,i);
	
	//3D, need to sort neighbors appropriately.
	//try LR-DU-BF, where F=front (high z) and B=back (low z)
	if (verbose == 1)
	  {
	    printf("\tPerforming Cartesian neighbor sorting: e1=%i, e2=%i, nx=%e, ny=%e, nz=%e\n",el1->getId(),el2->getId(),nx,ny,nz);
	  }
	if ((fabs(ny)+fabs(nz))<(2*eps) && fabs(nx) > eps)
	  {
	    //it's an x-normal face.
	    if (nx > 0)
	      {
		//el1 is -x, el2 is +x. el2 is right of el1
		_neighbors[el1num*N_N+1] = el2num;
		if (el2num>=0){ if(el2->getPartition()==_myid){ _neighbors[el2num*N_N+0] = el1num;}}
	      }
	    else //nx < 0
	      {
		//el1 is +x, el2 is -x. el2 is left of el1
		_neighbors[el1num*N_N+0] = el2num;
		// do the reverse if el2 belongs to me (and it's not el1num, ie el2num >=0 )
		if (el2num>=0){ if(el2->getPartition()==_myid){_neighbors[el2num*N_N+1] = el1num;}}
	      }
	    if (verbose == 1) {printf("\tThe interface is x-normal\n");}
	  } //end case for x-normal face
	else if ((fabs(nx)+fabs(nz))<(2*eps) && fabs(ny) > eps)
	  {
	    //It's a y-normal face
	    if (ny > 0)
	      {
		//el1 is -y, el2 is +y. el2 is up of el1
		_neighbors[el1num*N_N+3] = el2num;
		if (el2num>=0){ if(el2->getPartition()==_myid){ _neighbors[el2num*N_N+2] = el1num;}}
	      }
	    else //ny < 0
	      {
		//el1 is +y, el2 is -y. el2 is down of el1
		_neighbors[el1num*N_N+2] = el2num;
		if (el2num>=0){ if(el2->getPartition()==_myid){ _neighbors[el2num*N_N+3] = el1num;}}
	      }
	    if (verbose == 1) {printf("\tThe interface is y-normal\n");}
	  } //end case for y-normal face
	else if ((fabs(nx)+fabs(ny))<(2*eps) && fabs(nz) > eps)
	  {
	    //it's a z-normal face
	    if (nz > 0)
	      {
		//el1 is -z, el2 is +z. el2 is front of el1
		_neighbors[el1num*N_N+5] = el2num;
		if (el2num>=0){ if(el2->getPartition()==_myid){ _neighbors[el2num*N_N+4] = el1num;}}
	      }
	    else //nz < 0
	      {
		//el1 is +z, el2 is -z. el2 is back of el1
		_neighbors[el1num*N_N+4] = el2num;
		if (el2num>=0){ if(el2->getPartition()==_myid){ _neighbors[el2num*N_N+5] = el1num;}}
	      }
	    if (verbose == 1) {printf("\tThe interface is z-normal\n");}
	  } //end case for z-normal face
	else
	  {
	    printf("Error in neighbor creation!! The mesh was claimed Cartesian, but the normals do not back up the story.\n");
	    exit(1);
	  }
      } //end case structure for 3D
#endif
    } // end if cartesian
  } // end loop on interfaces
  if (verbose == 1) {printf("Left the interface loop, preparing to delete nn\n");}
  // Free some stuff
  delete[] nn;
}

simpleInterface::simpleInterface(int physicalTag)
{
  /*!
    \brief Constructor for an interface.
    \param[in] physicalTag the physical tag of the interface read from the mesh
  */
  _elements[0] = NULL;
  _elements[1] = NULL;
  _physicalTag = physicalTag;
  _closureId[0]= -1;
  _closureId[1]= -1;
}


void simpleInterface::BuildInterfaces(simpleMesh &mesh, std::vector<simpleInterface> &interfaces, int typeInterface, int typeElement, int nsides, int PeriCo2D)
{
  /*!
    \brief Generate the interfaces.
    \param[in] mesh the mesh
    \param[out] interfaces the interface vector
    \param[in] typeInterface the associated key number referencing the type of interface
    \param[in] typeElement the associated key number referencing the type of element
    \param[in] nsides number of sides to an element
    \param[in] PeriCo2D whether or not to make the HiOCFD5 vortex mesh correction

    \section Description
    This function builds all the interfaces. For each interface, it
    will find both elements on the right and left. If it's on a
    boundary, it will use the BC information to find the correct
    neighboring element.
    
    For MPI processes, there a bunch of conditions to take care of so
    things might get a bit convoluted. Hang on to your hats.
  */

  std::map<std::vector<int>, simpleInterface> interfacesMap;
  //
  // Read the pre-existing interfaces in the mesh (my partition)
  //
  int verbose = 0;
  if (verbose == 1) {printf("Entered BuildInterfaces\n");}
  //  int order = 1;
  //the mesh.getElements command looks only for interfaces, using typeInterface variable
  //that must be properly specified in main based on element type.
  const std::vector<simpleElement> preElements = mesh.getElements(typeInterface);
  if (verbose == 1){ printf("preElements size %ld\n", preElements.size());}
  for (int i = 0; i < preElements.size(); i++) {
    const simpleElement &el = preElements[i];
    std::vector<int> LOCALnodes;
    if (verbose == 1){printf("element %i\n", i);}
    for (int k = 0; k < el.getNbNodes(); k++) {
      LOCALnodes.push_back(el.getNode(k));
      if (verbose == 1){
      printf("     node (+1) %i\n",el.getNode(k)+1);
      printf("nodes[%d] = %d\n", k, LOCALnodes[k]);
      }
    }
    
    std::sort(LOCALnodes.begin(), LOCALnodes.end());
    if (verbose == 1)
      {
	printf("Just conducted sort operation on node indexing:\n");
	for (int k = 0; k < el.getNbNodes(); k++)
	  {
	    printf("nodes[%d] = %d\n", k, LOCALnodes[k]); }
      }
    //now, take the set of integers in 'nodes' and pair them with
    //the single physical tag stored for the interface
    interfacesMap[LOCALnodes] = simpleInterface(el.getPhysicalTag());
    if (verbose == 1) {printf("physical tag associated with preElement is %d\n", el.getPhysicalTag());}
  }
  //I think that is just to get all the boundary interfaces. Now, we have to loop 
  //through elements to find all the interior interfaces
  //  printf("Map size is %i after the preElement analysis\n", interfacesMap.size());

  //
  // Now loop through the elements in the mesh and find all the other
  // interfaces (in my partition)
  //
  if (verbose == 1){printf("Finished with boundary preElements, now looking for interior interfaces\n");}
  const std::vector<simpleElement> &elements = mesh.getElements(typeElement); 
  //Time to build the basis associated with the element:
  const polynomialBasis &basis = *polynomialBases::find(typeElement);
  const std::vector<std::vector<int> > &closures = basis.closures;
  if (verbose == 1)
    {
      printf("Just imported closures:same for all elements\n");
      for (int c = 0; c < closures.size(); c++)
	{
#ifdef TWOD
	  printf("closure[%d] is {%d,%d}\n",c,closures[c][0],closures[c][1]);
#endif
#ifdef THREED
	  printf("closure %d:\n",c);
	  for (int j = 0; j < closures[c].size(); j++)
	    {
	      printf("closure[%d][%d] = %d\n", c, j, closures[c][j]);
	    }
	  //printf("closure[%d] is {%d,%d,%d,%d}\n",c,closures[c][0],closures[c][1],closures[c][2],closures[c][3]);
#endif
	}
    }
  std::vector<int> nodes;
  for (int i = 0; i < elements.size(); i++) 
    {
      const simpleElement &el = elements[i];
      if (verbose == 1)
	{
	  printf("element %i\n", i);
	  printf("   id of the element is %i\n",(el).getId());
	  printf("   num of nodes is %i\n",(el).getNbNodes());
	  printf("   list of nodes:");
	  for(int p=0; p< el.getNbNodes();p++){
	    printf(" %i",el.getNode(p));
	  }
	  
	  printf("\n");
	  printf("The element's closure list, master nodes\n");
	}
      for (int j = 0; j < closures.size(); j++) 
      //for (int j = 0; j < basis.numFaces; j++) //I think I put this here for 3D
	{
	  nodes.clear();
	  //	  nodes.resize(0);
	  //the push_back command adds to size of nodes each time
	  // new value is necessary, so using the resize operation
	  // preemptively results in erroneous nodes values (a bunch of zeros)
	  //nodes.resize(pow(order+1,2)); //nodes per face
	  if (verbose == 1) {printf("j=%d:\n",j);}
	  const std::vector<int> &cl = closures[j];
	  for (int k = 0; k < cl.size(); k++) {
	    nodes.push_back(el.getNode(cl[k]));
	    if (verbose == 1) {printf("   node %i\n", nodes[k]);}
	  }
	  std::sort(nodes.begin(), nodes.end());
	  if (verbose == 1)
	    {
	      printf("Just conducted sort operation on node indexing:\n");
	      for (int k = 0; k < nodes.size(); k++)
		{
		  printf("nodes[%d] = %d\n", k, nodes[k]); }
	    }
	  //Okay: still inside the j loop, meaning we are working
	  //with the nodes of a specific face of a specific element

	  // If found in the map, use it. If not, add it to the map
	  simpleInterface &interfaceFound = interfacesMap[nodes];
	  if (interfaceFound._elements[0] == &el)
	    {
	      if (verbose == 1) {printf("This interface has already been claimed by this element, moving on\n");}
	    }
	  if (interfaceFound._elements[0] == NULL) 
	    {
	      if (verbose == 1) {printf("---element %i, (id=%d), j=%d: We've found a new interface, no known correspondence to node set---\n",i,(el).getId(),j);}
	      interfaceFound._elements[0] = &el;
	      interfaceFound._closureId[0] = j;
	      //farfield or reflective BC, copy my element to my neighbor
	      //if ((interfaceFound.getPhysicalTag()==2)||(interfaceFound.getPhysicalTag()==3)||(interfaceFound.getPhysicalTag()==4)){
	      //BR2 Edit: Add the 7 tag for the no-slip boundary to this case structure:
	      if ((interfaceFound.getPhysicalTag()==2)|| //This is zero-gradient boundary
		  (interfaceFound.getPhysicalTag()==3)|| //reflective tag
		  (interfaceFound.getPhysicalTag()==5)|| //was previously 4
		  (interfaceFound.getPhysicalTag()==9)|| //9 is the A inflow designation
		  (interfaceFound.getPhysicalTag()==11)|| //11 is the KJet inflow designation
		  (interfaceFound.getPhysicalTag()==13)|| //11 is the SubOut designation
		  (interfaceFound.getPhysicalTag()==101)|| //101 is the Homogeneous BC designation
		  (interfaceFound.getPhysicalTag()==7)) //no-slip tag
		{
		  if (verbose == 1) {printf("Tag indicates physical boundary, repeating el0 info for el1\n");}
		  interfaceFound._elements[1] = &el;
		  interfaceFound._closureId[1] = j;
		}
	      else{	  
		interfaceFound._elements[1] = NULL;
		interfaceFound._closureId[1] = -1;
		if (verbose == 1) {printf("\t\tSet interface's element1 and closure1 to null and -1\n");}
	      }
	//	printf("This has to be 8 times\n");
	    } //end case for _elemen[0]=NULL
	  else if (interfaceFound._elements[0] != &el) 
	    {
	      if (verbose == 1){printf("--The interface associated with these nodes has previously been claimed by a different element, specifically element id %d using closure id %d---\n",interfaceFound._elements[0]->getId(),interfaceFound._closureId[0]);}
	      if (interfaceFound._elements[1] == NULL) 
		{
		  
		  interfaceFound._elements[1] = &el;
		  interfaceFound._closureId[1] = j+nsides; 
		  if (verbose == 1) {printf("\t\tThere was room for a second element, so el %i (id=%d) is now _el[1] of this interface, using closure id %d\n", i,interfaceFound._elements[1]->getId(),interfaceFound._closureId[1]);}
		  //	  printf("This has to be 4 times\n");
		}
	      else if (interfaceFound._elements[1] != &el) {
		printf("error in interfaces creation when looking in my partition!!!\n");
		exit(1);
	      }
	    } //end it for _el[0] already claimed
	} //end j loop over the 2xfaces of element
    } //end loop over all the elements
  if (verbose > 0) {printf("\n\nConcluded the loop over all elements to find interior interfaces\n");}
#ifdef THREED
  {
    if (verbose == 1){printf("Progress report\n");}
    int CounterPhil = 0;
    for (std::map<std::vector<int>, simpleInterface>::iterator it = interfacesMap.begin(); it != interfacesMap.end(); ++it) {
      std::vector<int> nodes1 = it -> first;
      simpleInterface & interface1 = it -> second;
      if (verbose == 1)
	{
	  printf("Interface (%d) defined by nodes {%d,%d,%d,%d}, or GMSH: {%d,%d,%d,%d}\t\t", CounterPhil, nodes1[0], nodes1[1], nodes1[2],nodes1[3],nodes1[0]+1,nodes1[1]+1,nodes1[2]+1,nodes1[3]+1);
	}
      CounterPhil++;
      if (verbose == 1) {
	printf("Element 0 id is %d;\t",interface1._elements[0]->getId());
	if (interface1.getElement(1)==NULL)
	  {
	    printf("Element 1 remains unclaimed\n");
	  }
	else
	  {
	    printf("Element 1 id is %d\n",interface1._elements[1]->getId());
	  }
      }
    }
  }
#endif
  printf("\n");
  //
  // Now account for the periodic BC which are on my partition
  //
  const fullMatrix<double> &meshNodesCONST = mesh.getNodes();  
  //PEJ 11/05/2017: Mess with the node locations so I can run HiOCFD5 vortex problem
   fullMatrix<double> meshNodes; meshNodes.resize(meshNodesCONST.size1(), meshNodesCONST.size2());
  for (int i = 0; i < meshNodes.size1(); i++){
    for (int j = 0; j < meshNodes.size2(); j++){
      meshNodes(i,j) = meshNodesCONST(i,j); }}
  if (PeriCo2D > 0)
    {
      //HiOCFD5 vortex transport, where the periodic mesh edges do not
      //align properly. I must address this problem by directly altering
      //node coordinates.
      //Also, the integer PeriCorrection2D should be number of elements along wall
      std::vector<scalar> yLeft; yLeft.clear();
      std::vector<int> indexLeft; indexLeft.clear();
      std::vector<scalar> yRight; yRight.clear();
      std::vector<int> indexRight; indexRight.clear();
      std::vector<scalar> xBase; xBase.clear();
      std::vector<int> indexBase; indexBase.clear();
      std::vector<scalar> xTop; xTop.clear();
      std::vector<int> indexTop; indexTop.clear();
      /*
      scalar yLeft[PeriCo2D];
      scalar yRight[PeriCo2D];
      scalar xBase[PeriCo2D];
      scalar xTop[PeriCo2D];
      */
      scalar eps_mesh = pow(10,-12);
      scalar Lmesh = 0.1; //mesh-dependent, fix according to your desires
      scalar x0 = 0.0;
      scalar y0 = 0.0;
      for (int j = 0; j < meshNodes.size2(); j++)
	{
	  scalar xLocal = meshNodes(0,j);
	  scalar yLocal = meshNodes(1,j);
	  if (fabs(xLocal-x0)<eps_mesh)
	    {
	      yLeft.push_back(yLocal);
	      indexLeft.push_back(j);
	    }
	  if (fabs(yLocal-y0)<eps_mesh)
	    {
	      xBase.push_back(xLocal);
	      indexBase.push_back(j);
	    }
	  if (fabs(xLocal-(x0+Lmesh))<eps_mesh)
	    {
	      yRight.push_back(yLocal);
	      indexRight.push_back(j);
	    }
	  if (fabs(yLocal-(y0+Lmesh))<eps_mesh)
	    {
	      xTop.push_back(xLocal);
	      indexTop.push_back(j);
	    }
	}
      //Should now have WallRes count of x or y coordinate list along each wall.
      //Now, get the sorted order
      std::vector<int> orderLeft; orderLeft.resize(yLeft.size());
      std::vector<int> orderRight; orderRight.resize(yRight.size());
      std::vector<int> orderBase; orderBase.resize(xBase.size());
      std::vector<int> orderTop; orderTop.resize(xTop.size());
      for (int j1 = 0; j1 < orderLeft.size(); j1++) {
	  int index = 0;
	  for (int j2 = 0; j2 < orderLeft.size(); j2++)	{
	    if (yLeft[j2] < yLeft[j1])	{ index++; } }
	  orderLeft[j1] = index; }
      for (int j1 = 0; j1 < orderRight.size(); j1++) {
	  int index = 0;
	  for (int j2 = 0; j2 < orderRight.size(); j2++)	{
	    if (yRight[j2] < yRight[j1])	{ index++; } }
	  orderRight[j1] = index; }
      for (int j1 = 0; j1 < orderBase.size(); j1++) {
	  int index = 0;
	  for (int j2 = 0; j2 < orderBase.size(); j2++)	{
	    if (xBase[j2] < xBase[j1])	{ index++; } }
	  orderBase[j1] = index; }
      for (int j1 = 0; j1 < orderTop.size(); j1++) {
	  int index = 0;
	  for (int j2 = 0; j2 < orderTop.size(); j2++)	{
	    if (xTop[j2] < xTop[j1])	{ index++; } }
	  orderTop[j1] = index; }
      if (verbose > 0)
	{
	  //Let's have a look at node locations and the order arrays:
	  for (int j = 0; j < yLeft.size(); j++)
	    {
	      printf("yLeft[%d] = %f, order = %d, node index=%d\t||\t", j, yLeft[j], orderLeft[j], indexLeft[j]);
	      printf("yRight[%d] = %f, order = %d, node index = %d\n", j, yRight[j], orderRight[j], indexRight[j]);
	    }
	  for (int j = 0; j < xBase.size(); j++)
	    {
	      printf("xBase[%d] = %f, order = %d, node index = %d\t||\t", j, xBase[j], orderBase[j], indexBase[j]);
	      printf("xTop[%d] = %f, order = %d, node index = %d\n", j, xTop[j], orderTop[j], indexTop[j]);
	    }
	  
	}
      //Now, alter the right and top coordinates to eliminate hanging nodes
      for (int j1 = 0; j1 < orderRight.size(); j1++)
	{
	  //int place = orderRight[j];
	  for (int j2 = 0; j2 < orderRight.size(); j2++)
	    {
	      if (orderLeft[j2] == orderRight[j1])
		{
		  //We have found target y value
		  meshNodes(1, indexRight[j1]) = meshNodes(1, indexLeft[j2]);
		}
	    }
	}
      for (int j1 = 0; j1 < orderTop.size(); j1++)
	{
	  //int place = orderRight[j];
	  for (int j2 = 0; j2 < orderTop.size(); j2++)
	    {
	      if (orderBase[j2] == orderTop[j1])
		{
		  //We have found target x value
		  meshNodes(0, indexTop[j1]) = meshNodes(0, indexBase[j2]);
		}
	    }
	}
      if (verbose == 1)
    {
      printf("The nodes from getNodes, AFTER periodicity correction. nodes(direction, index):\n");
      for (int j = 0; j < meshNodes.size2(); j++)
	{
	  for (int i = 0; i < meshNodes.size1(); i++)
	    {
	      //these are physical node locations read in from the mesh file
	      printf("nodes(%d,%d)=%f,  ",i,j,meshNodes(i,j));} printf("\n"); }
    }
    }
  for (std::map<std::vector<int>, simpleInterface>::iterator it = interfacesMap.begin(); it != interfacesMap.end(); ++it) {
    std::vector<int> nodes1 = it -> first;
    simpleInterface & interface1 = it -> second;
#ifdef TWOD
    if (verbose == 1) {printf("The periodic BC loop: we are working on the interface defined by nodes {%d,%d}\t\t", nodes1[0], nodes1[1]);}
#endif
#ifdef THREED
    if (verbose == 1) {printf("The periodic BC loop: we are working on the interface defined by nodes {%d,%d,%d,%d}\t\t", nodes1[0], nodes1[1], nodes1[2], nodes1[3]);}
#endif
    if (verbose == 1)
      {
	printf("Element 0 id is %d;\t",interface1._elements[0]->getId());
	if (interface1.getElement(1)==NULL)
	  {
	    printf("Element 1 remains unclaimed\n");
	  }
	else
	  {
	    printf("Element 1 id is %d\n",interface1._elements[1]->getId());
	  }
      }
    // If it's a periodic BC and we haven't found a partner element yet
    if ((interface1.getPhysicalTag()==1)&&(interface1.getElement(1)==NULL)){
      if (verbose == 1) {printf("\t\tThe interface is tagged as a periodic face and el1 is unassigned\n");}
      // Now loop on all the other interfaces to find the match
      for (std::map<std::vector<int>, simpleInterface>::iterator it2 = interfacesMap.begin(); it2 != interfacesMap.end(); ++it2) {
  	std::vector<int> nodes2 = it2 -> first;
  	simpleInterface & interface2 = it2 -> second;

  	// If it's a periodic BC and we haven't found a partner element yet
  	if ((interface2.getPhysicalTag()==1)&&(interface2.getElement(1)==NULL)){
    //printf("nodes1={%d,%d,%d,%d}, nodes2={%d,%d,%d,%d}\n", nodes1[0],nodes1[1],nodes1[2],nodes1[3],nodes2[0],nodes2[1],nodes2[2],nodes2[3]);
  	  // Check if they are paired
  	  bool paired = pairPeriodic(meshNodes, nodes1, nodes2);
  	  if(paired){
  	    interface1._elements[1] = interface2.getElement(0);
  	    interface1._closureId[1] = interface2.getClosureId(0)+nsides;
  	    interface2._elements[1] = interface1.getElement(0);
  	    interface2._closureId[1] = interface1.getClosureId(0)+nsides;
	    if (verbose == 1) {printf("\t\tAchieved periodic pairing, partner element is el id %d\n",interface1._elements[1]->getId());}
	    //	    printf("\t
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
  if (verbose > 0){printf("\nInterface Pairing Check at end of BuildInterfaces\n");}
  for (std::map<std::vector<int>, simpleInterface>::iterator it = interfacesMap.begin(); it != interfacesMap.end(); ++it) {
    simpleInterface & interface = it -> second;
    if (verbose > 0)
      {
	printf("interface el1=%i\t\t",interface.getElement(0)->getId());
	printf("          el2=%i\n",interface.getElement(1)->getId());
      }
    if((interface.getElement(0)==NULL)||(interface.getElement(1)==NULL)){
      printf("error in interfaces creation !!!\n");
      printf("Failed on final check on element el=%i\n",interface.getElement(0)->getId());
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

void simpleInterface::FixClosures(simpleMesh &mesh, std::vector<simpleInterface> &interfaces, int typeInterface, int typeElement, int nsides, int PeriCo2D, scalar DomainLength)
{
  //PEJ 11/17/2017: Check that all interfaces are storing
  //proper closure information. If not, adjust the closure information
  printf("Entered FixClosures routine\n");
  int verbose = 0;
  const std::vector<simpleElement> &elements = mesh.getElements(typeElement); 
  const fullMatrix<double> &meshNodesCONST = mesh.getNodes();  
  //PEJ 11/05/2017: Mess with the node locations so I can run HiOCFD5 vortex problem
  fullMatrix<double> meshNodes; meshNodes.resize(meshNodesCONST.size1(), meshNodesCONST.size2());
  for (int i = 0; i < meshNodes.size1(); i++){
    for (int j = 0; j < meshNodes.size2(); j++){
      meshNodes(i,j) = meshNodesCONST(i,j); }}
  if (PeriCo2D > 0)
    {
      //HiOCFD5 vortex transport, where the periodic mesh edges do not
      //align properly. I must address this problem by directly altering
      //node coordinates.
      //Also, the integer PeriCorrection2D should be number of elements along wall
      std::vector<scalar> yLeft; yLeft.clear();
      std::vector<int> indexLeft; indexLeft.clear();
      std::vector<scalar> yRight; yRight.clear();
      std::vector<int> indexRight; indexRight.clear();
      std::vector<scalar> xBase; xBase.clear();
      std::vector<int> indexBase; indexBase.clear();
      std::vector<scalar> xTop; xTop.clear();
      std::vector<int> indexTop; indexTop.clear();
      /*
      scalar yLeft[PeriCo2D];
      scalar yRight[PeriCo2D];
      scalar xBase[PeriCo2D];
      scalar xTop[PeriCo2D];
      */
      scalar eps_mesh = pow(10,-12);
      scalar Lmesh = 0.1; //mesh-dependent, fix according to your desires
      scalar x0 = 0.0;
      scalar y0 = 0.0;
      for (int j = 0; j < meshNodes.size2(); j++)
	{
	  scalar xLocal = meshNodes(0,j);
	  scalar yLocal = meshNodes(1,j);
	  if (fabs(xLocal-x0)<eps_mesh)
	    {
	      yLeft.push_back(yLocal);
	      indexLeft.push_back(j);
	    }
	  if (fabs(yLocal-y0)<eps_mesh)
	    {
	      xBase.push_back(xLocal);
	      indexBase.push_back(j);
	    }
	  if (fabs(xLocal-(x0+Lmesh))<eps_mesh)
	    {
	      yRight.push_back(yLocal);
	      indexRight.push_back(j);
	    }
	  if (fabs(yLocal-(y0+Lmesh))<eps_mesh)
	    {
	      xTop.push_back(xLocal);
	      indexTop.push_back(j);
	    }
	}
      //Should now have WallRes count of x or y coordinate list along each wall.
      //Now, get the sorted order
      std::vector<int> orderLeft; orderLeft.resize(yLeft.size());
      std::vector<int> orderRight; orderRight.resize(yRight.size());
      std::vector<int> orderBase; orderBase.resize(xBase.size());
      std::vector<int> orderTop; orderTop.resize(xTop.size());
      for (int j1 = 0; j1 < orderLeft.size(); j1++) {
	  int index = 0;
	  for (int j2 = 0; j2 < orderLeft.size(); j2++)	{
	    if (yLeft[j2] < yLeft[j1])	{ index++; } }
	  orderLeft[j1] = index; }
      for (int j1 = 0; j1 < orderRight.size(); j1++) {
	  int index = 0;
	  for (int j2 = 0; j2 < orderRight.size(); j2++)	{
	    if (yRight[j2] < yRight[j1])	{ index++; } }
	  orderRight[j1] = index; }
      for (int j1 = 0; j1 < orderBase.size(); j1++) {
	  int index = 0;
	  for (int j2 = 0; j2 < orderBase.size(); j2++)	{
	    if (xBase[j2] < xBase[j1])	{ index++; } }
	  orderBase[j1] = index; }
      for (int j1 = 0; j1 < orderTop.size(); j1++) {
	  int index = 0;
	  for (int j2 = 0; j2 < orderTop.size(); j2++)	{
	    if (xTop[j2] < xTop[j1])	{ index++; } }
	  orderTop[j1] = index; }
      if (verbose > 0)
	{
	  //Let's have a look at node locations and the order arrays:
	  for (int j = 0; j < yLeft.size(); j++)
	    {
	      printf("yLeft[%d] = %f, order = %d, node index=%d\t||\t", j, yLeft[j], orderLeft[j], indexLeft[j]);
	      printf("yRight[%d] = %f, order = %d, node index = %d\n", j, yRight[j], orderRight[j], indexRight[j]);
	    }
	  for (int j = 0; j < xBase.size(); j++)
	    {
	      printf("xBase[%d] = %f, order = %d, node index = %d\t||\t", j, xBase[j], orderBase[j], indexBase[j]);
	      printf("xTop[%d] = %f, order = %d, node index = %d\n", j, xTop[j], orderTop[j], indexTop[j]);
	    }
	  
	}
      //Now, alter the right and top coordinates to eliminate hanging nodes
      for (int j1 = 0; j1 < orderRight.size(); j1++)
	{
	  //int place = orderRight[j];
	  for (int j2 = 0; j2 < orderRight.size(); j2++)
	    {
	      if (orderLeft[j2] == orderRight[j1])
		{
		  //We have found target y value
		  meshNodes(1, indexRight[j1]) = meshNodes(1, indexLeft[j2]);
		}
	    }
	}
      for (int j1 = 0; j1 < orderTop.size(); j1++)
	{
	  //int place = orderRight[j];
	  for (int j2 = 0; j2 < orderTop.size(); j2++)
	    {
	      if (orderBase[j2] == orderTop[j1])
		{
		  //We have found target x value
		  meshNodes(0, indexTop[j1]) = meshNodes(0, indexBase[j2]);
		}
	    }
	}
      if (verbose == 1)
    {
      printf("The nodes from getNodes, AFTER periodicity correction. nodes(direction, index):\n");
      for (int j = 0; j < meshNodes.size2(); j++)
	{
	  for (int i = 0; i < meshNodes.size1(); i++)
	    {
	      //these are physical node locations read in from the mesh file
	      printf("nodes(%d,%d)=%f,  ",i,j,meshNodes(i,j));} printf("\n"); }
    }
    }
  const polynomialBasis &basis = *polynomialBases::find(typeElement);
  const std::vector<std::vector<int> > &closures = basis.closures;
  for (int t = 0; t < interfaces.size(); t++)
    {
      if (verbose > 0) {printf("--Treating interface %d:--\n",t);}
      //Get the two elements
      const simpleElement* el0 = interfaces[t]._elements[0];
      const simpleElement* el1 = interfaces[t]._elements[1]; //may crash code when boundaries are involved
      //Get the closure ids of the two elements
      int j0 = interfaces[t]._closureId[0];
      int j1 = interfaces[t]._closureId[1];
      if (verbose > 0) {printf("\tomA=%d, omB=%d, closureA=%d, closureB=%d\n", el0->getId(), el1->getId(), j0,j1);}
      //Now, check if the global node indices corresponding to the two closures match
      std::vector<int> nodes0; nodes0.clear();
      std::vector<int> nodes1; nodes1.clear();
      std::vector<int> cl0 = closures[j0];
      std::vector<int> cl1 = closures[j1];
      for (int i = 0; i < cl0.size(); i++)
	{
	  nodes0.push_back(el0->getNode(cl0[i]));
	}
      for (int i = 0; i < cl1.size(); i++)
	{
	  nodes1.push_back(el1->getNode(cl1[i]));
	}
      //nodes0 and nodes1 now hold global node indices corresponding to the face.
      int match = 1; //assume match
      for (int i = 0; i < cl0.size(); i++)
	{
	  if (verbose > 0)
	    {
	      printf("\tnodes0[%d]=%d, nodes1[%d]=%d\n",i,nodes0[i],i,nodes1[i]);
	    }
	  if (nodes0[i] != nodes1[i])
	    {
	      if (verbose == 1)
		{
		  printf("\t\tNODES DO NOT MATCH; Correction possibly necessary depending on periodicity\n");
		}
	      //	    i = cl0.size(); //exit loop
	      match = 0;
	    }
	}
      //if match == 0, that means the node-closure mapping is wrong and we need
      //to use a different closure.
      //Assume el0 has proper closure and change e11 closure.
      //Also, if the edge is a periodic edge, don't mess with closure info just yet.
      //So, must determine if edge is periodic, in this case, the nodes indices
      //won't just be jumbled, but instead are different sets, and I will
      //need to check physical coordinates directly. In fact, I think that's what
      //I will do regardless
      //Only the first nodes's (nodes0[0] and nodes1[0]) physical location needs to be checked;
      //if there is a mismatch, then I flip the closure
      scalar x0 = meshNodes(0, nodes0[0]); scalar x1 = meshNodes(0, nodes1[0]);
      scalar y0 = meshNodes(1, nodes0[0]); scalar y1 = meshNodes(1, nodes1[0]);
      scalar z0 = meshNodes(2, nodes0[0]); scalar z1 = meshNodes(2, nodes1[0]);
      if (match == 0)
	{
	  //bool paired = pairPeriodic(meshNodes, nodes0, nodes1);
	  //if (paired){ 
	  // }
	  scalar x0 = meshNodes(0, nodes0[0]); scalar x1 = meshNodes(0, nodes1[0]);
	  scalar y0 = meshNodes(1, nodes0[0]); scalar y1 = meshNodes(1, nodes1[0]);
	  scalar z0 = meshNodes(2, nodes0[0]); scalar z1 = meshNodes(2, nodes1[0]);
	  /*if (paired)
	    {
	      PeriodicAdjustment(x0,y0,z0,x1,y1,z1);
	      }*/
	  scalar xdiff = fabs(x0-x1);
	  scalar ydiff = fabs(y0-y1);
	  scalar zdiff = fabs(z0-z1);
	  if (verbose > 0)
	    {
	      printf("\tentered match==0 loop: X0=(%f,%f,%f) and X1=(%f,%f,%f) and DIFF=(%f,%f,%f) and DomaiLength=%f\n",x0,y0,z0,x1,y1,z1,xdiff,ydiff,zdiff,DomainLength);
	    }
	  scalar tol = pow(10,-10);
	  //Treat for periodicity
	  if (fabs(xdiff-DomainLength) < tol)
	    {
	      //x-periodic edge
	      xdiff = xdiff - DomainLength;
	      if (verbose > 0)
		{
		  printf("\tx-periodic, new xdiff is %f\n",xdiff);
		}
	    }
	   if (fabs(ydiff-DomainLength) < tol)
	    {
	      //y-periodic edge
	      ydiff = ydiff - DomainLength;
	      if (verbose > 0)
		{
		  printf("\ty-periodic, new ydiff is %f\n",ydiff);
		}
	    }
	    if (fabs(zdiff-DomainLength) < tol)
	    {
	      //z-periodic edge
	      zdiff = zdiff - DomainLength;
	      if (verbose > 0)
		{
		  printf("\tz-periodic, new zdiff is %f\n",zdiff);
		}
	    }
	  if (xdiff + ydiff + zdiff > tol)
	    {
	      if (verbose > 0)
		{
		  printf("\tCurrent DIFF=(%f,%f,%f), so we are using the wrong closure, which is presently %d\n",xdiff,ydiff,zdiff,interfaces[t]._closureId[1]);
		}
	      //true mismatch between nodes; flip the closure for second element
	      interfaces[t]._closureId[1] = (interfaces[t]._closureId[1] + nsides) % (2*nsides);
	      if (verbose > 0)
		{
		  printf("\t t=%d:Altered el1 closure to %d\n",t,interfaces[t]._closureId[1]);
		}
	    }
	}
      /*
      bool paired = pairPeriodic(meshNodes, nodes0, nodes1);
      if (paired){ match = 1; }
      if (match == 0)
	{
	  interfaces[t]._closureId[1] = (interfaces[t]._closureId[1] + nsides) % (2*nsides);
	  if (verbose > 0)
	    {
	      printf("\t t=%d:Altered el1 closure to %d\n",t,interfaces[t]._closureId[1]);
	    }
	}
      */
    } //end loop over all interfaces
} //End subroutine


void simpleMesh::buildElementMap(int elem_type){
  /*!
    \brief Build a map relating element id to element index
    \param[in] elem_type the associated key number referencing that element
    \section Description
    Map for element ID -> element index in order of the U matrix
  */
  int verbose = 0;
  const std::vector<simpleElement> &elements = getElements(elem_type);
  for(int e = 0; e <elements.size(); e++){
    const simpleElement &el = elements[e];
    if (verbose == 1)
      {
	printf("e:%i, id:%i\n", e, el.getId());
      }
    _elementMap[el.getId()] = e;
  }
}
				      
void simpleMesh::buildCommunicators(int elem_type){

  /*!
    \brief Build the necessary maps and matrices to communicate between processes.
    \param[in] elem_type the associated key number referencing that element

    \section Description
    Create ghostElementMap. If el2 of an interface belongs to another
    partition, store a unique number. This will be used to access
    back columns of U/A in the limiting procedure (using the
    neighbors vector)
    
    Create ghostElementSend matrix containing the element index to
    send to other partitions, the number of that partition to send
    to, and its global id (as tag for communication).
    
    Create ghostElementRecv matrix containing the element index to
    store an element from another partition, the source partition of
    that element, and its global id (as tag for communication).
  */

  // Count the number of interfaces which have el1 and el2 on
  // different partitions and fill the ghostElementMap
  /*
  //Marc's original approach
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
      printf("myid=%d: just set _ghostElementMap[%d] = %d\n",_myid,el2->getId(),N_E + count);
      sharedCount[el2->getPartition()] = sharedCount[el2->getPartition()]++;
      count++;
    }
  }
  _N_ghosts = count;
  
  // initialize the vectors
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

      _ghostElementSend[count*3+0] = _elementMap.at(id1);
      _ghostElementSend[count*3+1] = partition2;
      _ghostElementSend[count*3+2] = id1;

      _ghostElementRecv[count*3+0] = _ghostElementMap.at(id2);
      _ghostElementRecv[count*3+1] = partition2;
      _ghostElementRecv[count*3+2] = id2;      
      
      printf("myid=%d: just set _ghostElemenSend[%d]=%d, _ghostElemenRecv[%d]=%d\n",_myid, count*3+0, _ghostElementSend[count*3+0] , count*3+0, _ghostElementRecv[count*3+0]);

      count++;	
    } // end if condition on partitions
  } // end for loop
  delete[] sharedCount;
  */
  //Phil's new approach, necessary to handle corners
  //in the partition boundaries
  int verbose = 0;
  int count = 0;
  const std::vector<simpleElement> &elements = getElements(elem_type);
  std::vector<int> dummySEND; dummySEND.clear();
  std::vector<int> dummyRECV; dummyRECV.clear();
  _ghostElement_Key.clear();
  _ghostElement_Val.clear();
  _ghostElement_HNe.clear(); //home-partition neighbor elId
  int N_E = elements.size(); // number of elements in my partition
  for(int i = 0; i < _interfaces.size(); i++) //loop over the partition's interfaces
    {
      const simpleInterface &face = _interfaces[i];  // get the interface
      const simpleElement *el1 = face.getElement(0); // get the element on one side
      const simpleElement *el2 = face.getElement(1); // get the element on the other side
      if( _myid != el2->getPartition()) //el2 is off-partition, time to jump through hoops
	{
	  int partition1 = el1->getPartition();
	  int id1 = el1->getId();
	  int partition2 = el2->getPartition();
	  int id2 = el2->getId();
	  
	  _ghostElement_Key.push_back(el2->getId()); //gmsh index of off-partition element
	  _ghostElement_Val.push_back(N_E+count); //corresponding ghost index on myid partition
	  _ghostElement_HNe.push_back(id1); //gmsh id of the on-partion neighbor of the ghost element.
	  if (verbose>0){	  printf("myid=%d: just set _ghostElement_HNe = %d, _ghostElement_Key = %d, _ghostElement_Val = %d\n",_myid, _ghostElement_HNe[count], _ghostElement_Key[count], _ghostElement_Val[count]);}
	  

	  //Corresponding Send and Receive array entries:
	  //Address of the flesh element sitting on myid partition:
	  dummySEND.push_back(_elementMap.at(id1));
	  //The opposing partition index:
	  dummySEND.push_back(partition2);
	  //The gmsh index of the flesh element sitting on myid partiion:
	  dummySEND.push_back(id1);
	  
	  //Now, that flesh element on myid processor acts
	  //as a ghost element for partition2
	  
	  //The partition2 address of the ghost element (lies between N_E and N_ghosts)
	  dummyRECV.push_back(_ghostElement_Val[count]);
	  //The opposing parition index:
	  dummyRECV.push_back(partition2);
	  //The gmsh index of the flesh element sitting on partition2
	  dummyRECV.push_back(id2);  
	  count++;
	} //close "if cross-partition" structure
    } //close i loop over interfaces

   _N_ghosts = count;
  
  // initialize the vectors
  _ghostElementSend = new int[3*_N_ghosts];
  _ghostElementRecv = new int[3*_N_ghosts];

  //copy the dummySEND and dummyRECV vectors in to ghostElementSend and ghostElementRecv
  for (int j = 0; j < 3*_N_ghosts; j++)
    {
      _ghostElementSend[j] = dummySEND[j];
      _ghostElementRecv[j] = dummyRECV[j];
    }
  if (verbose > 0)
    {
      for (int i = 0; i < _N_ghosts; i++)
	{
	  printf("myid=%d: _ghostElemenSend[%d]=%d, _ghostElemenRecv[%d]=%d\n",_myid, i*3+0, _ghostElementSend[i*3+0] , i*3+0, _ghostElementRecv[i*3+0]);
	}
    }
}
				      

// 
fullMatrix<scalar> simpleMesh::getElementCentroids(const int N_E, const int ncorners, const fullMatrix<scalar> XYZNodes){
  /*!
    \brief Calculate the element centroid coordinates
    \param[in] N_E number of elements
    \param[in] ncorners number of vertices
    \param[in] XYZNodes element nodal coordinates.
    \return Return a matrix which is N_E x D containing the element centroids
  */

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
  /*!
    \brief Figure out if the mesh is cartesian
    \param[in] typeElement the type of element we are dealing with
    \param[in] elem_type the associated key number referencing that element
    \return true if it's a cartesian mesh
  */
  
  
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
  else if (typeElement=="hex"){
    //PEJ addition 09/15/2017
    const std::vector<simpleElement> &elements = getElements(elem_type);
    const simpleElement &el = elements[0];
    //Get the eight corner nodes of the first element in the mesh
    scalar eps = pow(10,-9);
    scalar x[8];
    scalar y[8];
    scalar z[8];
    for (int j = 0; j < 8; j++)
      {
	x[j] = (scalar)_nodes(0,el.getNode(j));
	y[j] = (scalar)_nodes(1,el.getNode(j));
	z[j] = (scalar)_nodes(2,el.getNode(j));
      }
    //What defines a Cartesian element?
    //4 nodes with identical x
    //4 nodes with identical y
    //4 nodes with identical z
    int xcount = 0;
    int ycount = 0;
    int zcount = 0;
    scalar xref = x[0];
    scalar yref = y[0];
    scalar zref = z[0];
    for (int j = 0; j < 8; j++)
      {
	if (fabs(xref - x[j]) < eps)
	  {
	    xcount++;
	  }
	if (fabs(yref - y[j]) < eps)
	  {
	    ycount++;
	  }
	if (fabs(zref - z[j]) < eps)
	  {
	    zcount++;
	  }
      }
    printf("End of check for 3D Cartesian mesh: xcount=%d, ycount=%d, zcount=%d\n",xcount,ycount,zcount);
    if (xcount == 4 && ycount == 4 && zcount == 4)
      {
	cartesian = true;
	printf("IDENTIFIED CARTESIAN MESH\n");
      }
  }

  // set the local bool
  _cartesian = cartesian;
  
  return cartesian;
}


void simpleMesh::setDx(const int N_N, const int N_E, const fullMatrix<scalar> &XYZCen, const fullMatrix<scalar> &XYZNodes){
  /*!
    \brief Sets the minimum delta x in the mesh
    \param[in] N_N number of neighbor of one element
    \param[in] N_E number of elements
    \param[in] XYZCen element centroid coordinates
    \param[in] XYZNodes element nodal coordinates.
    
    \section Description
    This function finds the minium Dx of all the elements in the
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
  
#elif THREED
  // distance from point to surface?
  //Sorry Marc, I don't know either. To start, it's going to be
  //the minimum edge length.
  //Future reference: there is a point to surface distance routine
  //in recoverytools, could bring it in later
  //if I want to be more correct on deformed elements
  if (N_N != 6)
    {
      printf("CATASTROPHE!!! getDx function not ready for non-hex in 3D\n");
    }
  scalar dmin;
  scalar dx;
  scalar dy;
  scalar dz;
  //Initialize dmin as distance between first two nodes of first element
  dx = XYZNodes(0, 0*D+0) - XYZNodes(1, 0*D+0);
  dy = XYZNodes(0, 0*D+1) - XYZNodes(1, 0*D+1);
  dz = XYZNodes(0, 0*D+2) - XYZNodes(1, 0*D+2);
  dmin = sqrt(dx*dx + dy*dy + dz*dz);
  for (int om = 0; om < N_E; om++)
    {
      for (int j1 = 0; j1 < 8; j1++) //8 for the # of corner nodes per element
	{
	  for (int j2 = 0; j2 < 8; j2++)
	    {
	      if (j1 != j2)
		{
		  dx = XYZNodes(j1, om*D+0) - XYZNodes(j2, om*D+0);
		  dy = XYZNodes(j1, om*D+1) - XYZNodes(j2, om*D+1);
		  dz = XYZNodes(j1, om*D+2) - XYZNodes(j2, om*D+2);
		  scalar distLocal = sqrt(dx*dx + dy*dy + dz*dz);
		  dmin = fmin(dmin, distLocal);
		}
	    }
	}
    }
    printf("Exiting the 3D getDx routine: dmin = %f\n", dmin);
  _Dx = dmin;
#endif
    

}
inline scalar simpleMesh::distance_to_edge(scalar x0,scalar y0,scalar x1,scalar y1,scalar x2,scalar y2){
  /*!
    \brief Distance from center to an edge
    \param x0 x-coordinate of cell center
    \param y0 y-coordinate of cell center
    \param x1 x-coordinate of first node defining the edge
    \param y1 y-coordinate of first node defining the edge
    \param x2 x-coordinate of secon node defining the edge
    \param y2 y-coordinate of second node defining the edge
    \return distance
    \section Description
    from http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
  */
  return fabs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}


