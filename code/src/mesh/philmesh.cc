#include "philmesh.h"

static fullMatrix<double> gmshGeneratePointsTriangle(int order, bool serendip)
{
  int nbPoints = serendip ? 3 * order : (order + 1) * (order + 2) / 2;
  fullMatrix<double> point(nbPoints, 2);

  point(0, 0) = 0;
  point(0, 1) = 0;

  if (order > 0) {
    double dd = 1. / order;

    point(1, 0) = 1;
    point(1, 1) = 0;
    point(2, 0) = 0;
    point(2, 1) = 1;

    int index = 3;

    if (order > 1) {

      double ksi = 0;
      double eta = 0;

      for (int i = 0; i < order - 1; i++, index++) {
        ksi += dd;
        point(index, 0) = ksi;
        point(index, 1) = eta;
      }

      ksi = 1.;

      for (int i = 0; i < order - 1; i++, index++) {
        ksi -= dd;
        eta += dd;
        point(index, 0) = ksi;
        point(index, 1) = eta;
      }

      eta = 1.;
      ksi = 0.;

      for (int i = 0; i < order - 1; i++, index++) {
        eta -= dd;
        point(index, 0) = ksi;
        point(index, 1) = eta;
      }

      if (order > 2 && !serendip) {
        fullMatrix<double> inner = gmshGeneratePointsTriangle(order - 3, serendip);
        inner.scale(1. - 3. * dd);
        inner.add(dd);
        point.copy(inner, 0, nbPoints - index, 0, 2, index, 0);
      }
    }
  }
  return point;
}

static fullMatrix<double> gmshGeneratePointsQuad(int order, bool serendip)
{
  int nbPoints = serendip ? order*4 : (order+1)*(order+1);
  fullMatrix<double> point(nbPoints, 2);

  if (order > 0) {
    point(0, 0) = -1;
    point(0, 1) = -1;
    point(1, 0) = 1;
    point(1, 1) = -1;
    point(2, 0) = 1;
    point(2, 1) = 1;
    point(3, 0) = -1;
    point(3, 1) = 1;

    if (order > 1) {
      int index = 4;
      const static int edges[4][2]={{0,1},{1,2},{2,3},{3,0}};
      for (int iedge=0; iedge<4; iedge++) {
        int p0 = edges[iedge][0];
        int p1 = edges[iedge][1];
        for (int i = 1; i < order; i++, index++) {
          point(index, 0) = point(p0, 0) + i*(point(p1,0)-point(p0,0))/order;
          point(index, 1) = point(p0, 1) + i*(point(p1,1)-point(p0,1))/order;
        }
      }
      if (order > 2 && !serendip) {
        fullMatrix<double> inner = gmshGeneratePointsQuad(order - 2, false);
        inner.scale(1. - 2./order);
        point.copy(inner, 0, nbPoints - index, 0, 2, index, 0);
      }
    }
  }
  else {
    point(0, 0) = 0;
    point(0, 1) = 0;
  }
  return point;
}

fullMatrix<double> PeriFixNodes(int PeriCo2D, const fullMatrix<double> &nodesCONST)
{
  /*!
    \brief function to alter boundary nodes in HiOCFD5 vortex meshes
    to achieve conformality along periodic boundary interfaces
    \param[in] PeriCo2D whether or not to execure periodicity fix
    \param[in] nodesConst the stored mesh nodes which I have unfortunately been unable to modify
   */
  int verbose = 0;
  if (verbose > 0) {printf("Entered PeriFixNodes\n");}
  //Step 1: Size the output array and set it equal to the input array
  fullMatrix<double> nodes; //the output
  nodes.resize(nodesCONST.size1(), nodesCONST.size2());
  for (int i = 0; i < nodes.size1(); i++){
    for (int j = 0; j < nodes.size2(); j++){
      nodes(i,j) = nodesCONST(i,j); }}
  
  //Step 2: If signalled, then execute the periodicity fix on the nodes.
  //If not signalled, the nodes output is same as input.
  if (PeriCo2D > 0)
    {
      if (verbose > 0)
	{
	  printf("PeriCo2D>0 in PeriFixNodes, se we are executing periodicity fix\n");
	}
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
      scalar x0 = 0.0; //mesh-dependent
      scalar y0 = 0.0; //meh-dependent
      for (int j = 0; j < nodes.size2(); j++)
	{
	  scalar xLocal = nodes(0,j);
	  scalar yLocal = nodes(1,j);
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
		  nodes(1, indexRight[j1]) = nodes(1, indexLeft[j2]);
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
		  nodes(0, indexTop[j1]) = nodes(0, indexBase[j2]);
		}
	    }
	}
      //New issue: for p>1, I need to adjust the nodes of the elements
      //on right and top boundaries, can't just change the rightmost or topmost node.
      
      if (verbose == 1)
	{
	  printf("The nodes from getNodes, AFTER periodicity correction. nodes(direction, index):\n");
	  for (int j = 0; j < nodes.size2(); j++)
	    {
	      for (int i = 0; i < nodes.size1(); i++)
		{
		  //these are physical node locations read in from the mesh file
		  printf("nodes(%d,%d)=%f,  ",i,j,nodes(i,j));} printf("\n"); }
	}
    }
  return nodes;
}



fullMatrix<double> PeriFixNodes_HiOMesh(int PeriCo2D, fullMatrix<double> nodes, /*const std::vector<simpleElement> &elements*/ std::vector<std::vector<int> > elemNodes, int elem_type, int order, int N_N)
{
  /*!
    \brief 2nd function to mesh with node locations; this one corrects the non-vertex nodes in
    \p>1 elements
   */
  /*
    I don't think I can pass &elements in to this array,
    so I'm settling just for the node-node mapping
   */
  //Step 0: set output equal to input nodes
  int verbose = 0;
  fullMatrix<double> NewNodes;
  NewNodes.resize(nodes.size1(), nodes.size2());
  for (int i = 0; i < nodes.size1(); i++){
    for (int j = 0; j < nodes.size2(); j++){
      NewNodes(i,j) = nodes(i,j);
    }}

  if (PeriCo2D > 0)
    {
      if (order < 2)
	{
	  //do nothing, vertex nodes are already in place
	}
      else 
	{
	  //Step 1: get the element addresses of all elements
	  //that border the right or top boundary.
	  scalar eps_mesh = pow(10,-12);
	  scalar Lmesh = 0.1; //mesh-dependent, fix according to your desires
	  scalar x0 = 0.0; //mesh-dependent
	  scalar y0 = 0.0; //meh-dependent
	  int N_E = elemNodes.size();
	  if (verbose > 0) {printf("int PeriFixNodes_HiOMesh: N_E=%d\n", N_E);}
	  std::vector<int> elemBound; //element addresses, top/right boundaries
	  elemBound.clear();
	  for (int om = 0; om < N_E; om++) //element index
	    {
	      int N_s =  elemNodes[om].size();
	      //need separate k loops for x/y because of how I break the k loops when boundary element is found
	      for (int k = 0; k < N_s; k++) //local node index
		{
		  int J = elemNodes[om][k]; //global node index
		  scalar xLocal = nodes(0,J);
		  scalar yLocal = nodes(1,J);
		  //Query the nodes array to check if this element is a right/top boundary element
		  if (fabs(xLocal - (x0+Lmesh)) < eps_mesh)
		    {
		      //It's a right-boundary element.
		      elemBound.push_back(om);
		      k = N_s;
		    }
		  else if (fabs(yLocal - (y0+Lmesh)) < eps_mesh)
		    {
		      //It's a top-boundary element.
		      elemBound.push_back(om);
		      k = N_s;
		    }
		}
	    }
	  
	  //Now, for each element on either the top or right boundary,
	  //distribute nodes linearly under assumption that the vertex nodes
	  //(first N_N nodes) are properly located already.
	  //Get some reference element information:
	  //const polynomialBasis *basis  = polynomialBases::find(elem_type);  // for the element
	  //const std::vector<std::vector<int> > &closures = basis->closures;
	  fullMatrix<double> RefCords;
	  
	  if (N_N == 4)
	    {
	      //0 means "not serindiputy element"
	      RefCords = gmshGeneratePointsQuad(order, 0);
	      int N_s = RefCords.size1(); 
	      if (verbose > 0) {
		printf("In PeriFix_HiOMesh: N_s = %d\n", N_s); 
		printf("The reference coordinates:\n");
		for (int k = 0; k < N_s; k++)
		  {
		    printf("node %d: xi=%f, eta=%f\n",k,RefCords(k,0),RefCords(k,1));
		  }
	      }
	      //Wherever a node is on the reference element,
	      //that location should be reflected in NewNodes.
	      for (int i = 0; i < elemBound.size(); i++)
		{
		  //Get the vertex nodes for the element:
		  //elemNodes[om] = collection of N_s global node indices corresponding to element "om"
		  scalar x0 = nodes(0, elemNodes[elemBound[i]] [0]); scalar y0 = nodes(1, elemNodes[elemBound[i]] [0]); 
		  scalar x1 = nodes(0, elemNodes[elemBound[i]] [1]); scalar y1 = nodes(1, elemNodes[elemBound[i]] [1]); 
		  scalar x2 = nodes(0, elemNodes[elemBound[i]] [2]); scalar y2 = nodes(1, elemNodes[elemBound[i]] [2]); 
		  scalar x3 = nodes(0, elemNodes[elemBound[i]] [3]); scalar y3 = nodes(1, elemNodes[elemBound[i]] [3]); 
		  if (verbose > 0) {printf("Altering nodes of element %d: vertices are (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n",elemBound[i], x0,y0, x1,y1, x2,y2, x3,y3);}
		  for (int k = N_N; k < N_s; k++) //only mess with interior nodes
		    {
		      //Get the x/y coordinate of the solution point based on xi/eta coordinate and physical vertex coordinates.
		      //xiBi and etaBi come from bi-unit square domain; my formula for xNew/yNew uses unit domain
		      scalar xiBi = RefCords(k,0); scalar etaBi = RefCords(k,1);
		      scalar xi = 0.5 * (xiBi+1.0); scalar eta = 0.5*(etaBi+1.0);
		      scalar xNew = 0.0; scalar yNew = 0.0;
		      if (N_N == 4)
			{
			  xNew += (1.0-xi)*(1.0-eta)*x0; yNew += (1.0-xi)*(1.0-eta)*y0; 
			  xNew += (xi)    *(1.0-eta)*x1; yNew += (xi)    *(1.0-eta)*y1; 
			  xNew += (xi)    *(eta)    *x2; yNew += (xi)    *(eta)    *y2; 
			  xNew += (1.0-xi)*(eta)    *x3; yNew += (1.0-xi)*(eta)    *y3; 
			}
		      //Communicate new x/y coordinates to the NewNodes array
		      NewNodes(0, elemNodes[elemBound[i]] [k]) = xNew; NewNodes(1, elemNodes[elemBound[i]] [k]) = yNew;
		      if (verbose > 0) {printf("\t\tAltered local node %d, master index %d: old (x,y) = (%f,%f), new (x,y) = (%f, %f)\n", k, elemNodes[elemBound[i]][k], nodes(0, elemNodes[elemBound[i]] [k]), nodes(1, elemNodes[elemBound[i]] [k]), NewNodes(0, elemNodes[elemBound[i]] [k]), NewNodes(1, elemNodes[elemBound[i]] [k]));}
		    }
		}
	    } //end quadrilateral (N_N==4) case
	  if (N_N == 3)
	    {
	      //0 means "not serindipity element"
	      RefCords = gmshGeneratePointsTriangle(order, 0);
	      int N_s = RefCords.size1(); 
	      if (verbose > 0) {printf("In PeriFix_HiOMesh: N_s = %d\n", N_s); }
	      //Wherever a node is on the reference element,
	      //that location should be reflected in NewNodes.
	      for (int i = 0; i < elemBound.size(); i++)
		{
		  //Get the vertex nodes for the element:
		  scalar x0 = nodes(0, elemNodes[elemBound[i]] [0]); scalar y0 = nodes(1, elemNodes[elemBound[i]] [0]); 
		  scalar x1 = nodes(0, elemNodes[elemBound[i]] [1]); scalar y1 = nodes(1, elemNodes[elemBound[i]] [1]); 
		  scalar x2 = nodes(0, elemNodes[elemBound[i]] [2]); scalar y2 = nodes(1, elemNodes[elemBound[i]] [2]); 
	      
		  for (int k = N_N; k < N_s; k++) //only mess with interior nodes
		    {
		      //Get the x/y coordinate of the solution point based on xi/eta coordinate and physical vertex coordinates
		      scalar xi = RefCords(k,0); scalar eta = RefCords(k,1);
		      scalar xNew = 0.0; scalar yNew = 0.0;
	
		      xNew += x0;            yNew += y0;
		      xNew += (x1-x0) * xi;  yNew += (y1-y0) * xi;
		      xNew += (x2-x0) * eta; yNew += (y2-y0) * eta;
		  
		      //Communicate new x/y coordinates to the NewNodes array
		      NewNodes(0, elemNodes[elemBound[i]] [k]) = xNew; NewNodes(1, elemNodes[elemBound[i]] [k]) = yNew;
		    }
		}
	    } //end simplex (N_N == 3) case
	} //else "p>1"
    } //if "periCo2D>1"

  return NewNodes;
}
