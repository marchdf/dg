// Mesh for shock-drop problem.
// Specify the number of cells per radius, the radius, the lengths of the domain.
// Generates a mesh for the whole domain
// Does some grid stretching in the y-direction
/*

  (1)------------------------------------(2)
   |                                      |
   |                                      |
  (8)------------------------------------(3)
   |                                      |
   |                                      |
  (7)------------------------------------(4)
   |                                      |
   |                                      |
  (6)------------------------------------(5)

 */



lc = 0.0003;   // longueur caractéristique

// Domain lengths
cells_per_radius = 125;
radius = 1.0;
LxL = 7*radius; // length to the left of the drop
LxR = 16*radius; // length to the right of the drop

LyA = 4*radius;  // total length in y with uniform spacing
LyB = 6*radius;  // top/bottom length in y with grid stretching

// Sizes
dx = radius/cells_per_radius;
Nx=(LxL+LxR)/dx+1;
dy=dx;
NyA=LyA/dy+1;

// Solve for NyB assume geometric progression of rate r
r = 1.05;
NyB = Floor(Log(1-(1-r)*LyB/dy)/Log(r)+1);
//NyB=LyB/dy+1;
Printf("Number of nodes in stretch zone = %f",NyB);
  
// Définition des Points
Point(1) = { -LxL,   LyA/2+LyB,0,lc};  // upper left
Point(2) = {  LxR,   LyA/2+LyB,0,lc};   // upper right
Point(3) = {  LxR,   LyA/2,0,lc}; // middle upper right
Point(4) = {  LxR,  -LyA/2,0,lc}; // middle lower right
Point(5) = {  LxR, -(LyA/2+LyB),0,lc};  // lower right
Point(6) = { -LxL, -(LyA/2+LyB),0,lc}; // lower left
Point(7) = { -LxL,  -LyA/2,0,lc}; // middle lower left
Point(8) = { -LxL,   LyA/2,0,lc}; // middle upper left

// Définition des Lignes
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};

Line(9)  = {3,8};
Line(10) = {4,7};

// Définition de la Surface
Line Loop(55) = {1,2,9,8};
Line Loop(66) = {-9,3,10,7};
Line Loop(77) = {-10,4,5,6};
Physical Line(2) = {1,2,3,4,5,6,7,8};

Plane Surface(7) = {55};
Plane Surface(8) = {66};
Plane Surface(9) = {77};
Physical Surface(10) = {7,8,9};

Transfinite Line{1,5,9,10}=Nx; // x-resolution
Transfinite Line{-2,8,4,-6}=NyB Using Progression r; // y-resolution
Transfinite Line{3,7}=NyA; // y-resolution

Transfinite Surface(7)={1,2,3,8};
Transfinite Surface(8)={8,3,4,7};
Transfinite Surface(9)={7,4,5,6};

Mesh.RecombineAll=1;

Mesh.Partitioner = 1; // 1 = Chaco, 2 = Metis

// Chaco
Mesh.ChacoGlobalMethod = 1;
Mesh.ChacoArchitecture = 2;

//Metis
Mesh.MetisAlgorithm = 3;
Mesh.MetisEdgeMatching = 2;

// Partition with plugin
Plugin(SimplePartition).NumSlices = 4;
Plugin(SimplePartition).Run;
