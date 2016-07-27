// 
/*

  (1)-------(2)
   |         | 
   |         | 
   |         | 
   |         | 
   |         | 
   |         | 
  (4)-------(3)

 */

// Characterisic length
lc = 0.3;  

// Domain lengths
Lx = 1;
LyA = 2*Lx; // above the interface
LyB = 2*Lx; // below the interface

Nx=20+1;
dx=Lx/(Nx-1);
dy=dx;
Ny=(LyA+LyB)/dy+1;

// Définition des Points
Point(1) = { 0,  LyA,0,lc}; 
Point(2) = { Lx, LyA, 0,lc};  
Point(3) = { Lx,-LyB, 0,lc};
Point(4) = { 0, -LyB, 0,lc};

// Définition des Lignes
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Surface definitions
Line Loop(5) = {1,2,3,4};

// Physicals on the line (BC)
Physical Line(1) = {2,4};
Physical Line(2) = {1,3};

// Create the surfaces
Plane Surface(7) = {5};
Physical Surface(88) = {7};

// mesh part
Transfinite Line{1,3}=Nx; // x-resolution
Transfinite Line{2,4}=Ny; // y-resolution

Transfinite Surface(7)={1,2,3,4};

Mesh.RecombineAll=1;

Mesh.Partitioner = 1; // 1 = Chaco, 2 = Metis

// Chaco
Mesh.ChacoGlobalMethod = 1;
Mesh.ChacoArchitecture = 2;

//Metis
Mesh.MetisAlgorithm = 3;
Mesh.MetisEdgeMatching = 2;
