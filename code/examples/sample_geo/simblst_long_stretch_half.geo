// Very long version of simblst.geo with grid stretching, half the domain width

/*

  (1)-------(2)
   |         | 
   |         | 
   |         | 
  (8)-------(3) 
   |         | 
   |         |
   |         |
   |         |
   |         |
   |         |
   |         |
   |         | 
  (7)-------(4)
   |         | 
   |         | 
   |         | 
  (6)-------(5)

 */

// Characterisic length
lc = 0.3;  

// Domain lengths
lambda = 1;
Lx = 0.5;
LyA = 20*lambda; // above the interface
LyB = 30*lambda; // below the interface

Nx=20+1;
dx=Lx/(Nx-1);
dy=dx;
Ny=(LyA+LyB)/dy+1;

// stretching above and below the main domain
LS = 10*lambda; 

// Solve for NyS assume geometric progression of rate r
rS = 1.08;
NyS = Floor(Log(1-(1-rS)*LS/dy)/Log(rS)+1);
Printf("Number of nodes in stretch zone (y) = %f",NyS);  

// Définition des Points
Point(1) = { 0,  LyA+LS,0,lc}; 
Point(2) = { Lx, LyA+LS,0,lc}; 
Point(3) = { Lx, LyA, 0,lc};  
Point(4) = { Lx,-LyB, 0,lc};
Point(5) = { Lx,-LyB-LS, 0,lc};
Point(6) = { 0,-LyB-LS, 0,lc};
Point(7) = { 0,-LyB, 0,lc};
Point(8) = { 0, LyA, 0,lc};

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

// Surface definitions
Line Loop(5) = {1,2,9,8};
Line Loop(6) = {-9,3,10,7};
Line Loop(7) = {-10,4,5,6};

// Physicals on the line (BC)
Physical Line(2) = {1,5};
Physical Line(3) = {2,3,4,6,7,8};

// Create the surfaces
Plane Surface(7) = {5};
Plane Surface(8) = {6};
Plane Surface(9) = {7};
Physical Surface(88) = {7,8,9};

// mesh part
Transfinite Line{1,9,10,5}=Nx; // x-resolution
Transfinite Line{3,7}=Ny; // y-resolution

// stretching
Transfinite Line{-2,8,4,-6} = NyS Using Progression rS; 


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
