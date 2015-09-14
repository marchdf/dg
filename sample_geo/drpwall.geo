// Mesh for drop-wall problem.
// Specify the number of cells per radius, the radius, the lengths of the domain.
// Generates a mesh for the whole domain
// Does some grid stretching in the x-y directions (away from the wall)
/*

  (1)----------------(2)------------(3)
   |                  |              |
   |                  |              |
  (10)---------------(11)-----------(4)
   |                  |              |  wall
   |                  |              |
  (9)----------------(12)-----------(5)
   |                  |              |
   |                  |              |
  (8)----------------(7)------------(6)

 */



lc = 0.0003;   // longueur caractéristique

// Domain lengths
cells_per_radius = 75;
radius = 1.0;
LxW = 6*radius; // length in x for uniform spacing (near the wall)
LxF = 10*radius; // length in x for stretch spacing (far from wall)

LyA = 6*radius;  // total length in y with uniform spacing
LyB = 10*radius;  // top/bottom length in y with grid stretching

// Sizes
dx = radius/cells_per_radius;
NxW=(LxW)/dx+1;
dy=dx;
NyA=LyA/dy+1;

// Solve for NyB assume geometric progression of rate r
rB = 1.08;
NyB = Floor(Log(1-(1-rB)*LyB/dy)/Log(rB)+1);
Printf("Number of nodes in stretch zone (y) = %f",NyB);

// Solve for NxF assume geometric progression of rate r
rF = 1.08;
NxF = Floor(Log(1-(1-rF)*LxF/dx)/Log(rF)+1);
Printf("Number of nodes in stretch zone (x) = %f",NxF);

// Définition des Points
Point(1) = { -(LxF+LxW),   LyA/2+LyB,0,lc};
Point(2) = { -(LxW),       LyA/2+LyB,0,lc};
Point(3) = {  0,           LyA/2+LyB,0,lc}; 
Point(4) = {  0,           LyA/2,0,lc}; 
Point(5) = {  0,          -LyA/2,0,lc};
Point(6) = {  0,          -(LyA/2+LyB),0,lc};
Point(7) = { -(LxW),      -(LyA/2+LyB),0,lc}; 
Point(8) = { -(LxF+LxW),  -(LyA/2+LyB),0,lc};
Point(9) = { -(LxF+LxW),  -LyA/2,0,lc};
Point(10)= { -(LxF+LxW),   LyA/2,0,lc};
Point(11)= { -(LxW),       LyA/2,0,lc};
Point(12)= { -(LxW),      -LyA/2,0,lc};

// Définition des Lignes
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10)= {10,1};

Line(11) = {2,11};
Line(12) = {11,12};
Line(13) = {12,7};
Line(14) = {4,11};
Line(15) = {11,10};
Line(16) = {5,12};
Line(17) = {12,9};


// Définition de la Surface
Line Loop(44) = {1,11,15,10};
Line Loop(55) = {2,3,14,-11};
Line Loop(66) = {-15,12,17,9};
Line Loop(77) = {-14,4,16,-12};
Line Loop(88) = {-17,13,7,8};
Line Loop(99) = {-16,5,6,-13};
Physical Line(2) = {1,2,6,7,8,9,10};
Physical Line(3) = {3,4,5};

Plane Surface(4) = {44};
Plane Surface(5) = {55};
Plane Surface(6) = {66};
Plane Surface(7) = {77};
Plane Surface(8) = {88};
Plane Surface(9) = {99};
Physical Surface(10) = {4,5,6,7,8,9};

Transfinite Line{4,12,-9} = NyA;
Transfinite Line{2,-14,-16,-6} = NxW;
Transfinite Line{-1,15,17,7} = NxF Using Progression rF;
Transfinite Line{-3,-11,10,5,13,-8} = NyB Using Progression rB;

Transfinite Surface(4)={1,2,11,10};
Transfinite Surface(5)={2,3,4,11};
Transfinite Surface(6)={10,11,12,9};
Transfinite Surface(7)={11,4,5,12};
Transfinite Surface(8)={9,12,7,8};
Transfinite Surface(9)={12,5,6,7};

Mesh.RecombineAll=1;

Mesh.Partitioner = 1; // 1 = Chaco, 2 = Metis

// Chaco
Mesh.ChacoGlobalMethod = 1;
Mesh.ChacoArchitecture = 2;

//Metis
Mesh.MetisAlgorithm = 3;
Mesh.MetisEdgeMatching = 2;

// Partition with plugin
//Plugin(SimplePartition).NumSlices = 4;
//Plugin(SimplePartition).Run;
