//
// Long shock-tube problem (mostly used for simblst)
//
//


lc = 0.3;   // longueur caractéristique

// Domain lengths
Lx = 1;
LyA = 20*Lx; // above the interface
LyB = 30*Lx; // below the interface

Nx=50+1;
dx=Lx/(Nx-1);
dy=dx;
Ny=(LyA+LyB)/dy+1;

// Nx=256+1;
// dx=Lx/(Nx-1);
// dy=dx;
// Ny=(LyA+LyB)/dy+1;

// Définition des Points
Point(1) = { 0,  LyA,0,lc};  // upper left
Point(2) = { Lx, LyA,0,lc};   // upper right
Point(3) = { Lx,-LyB, 0,lc};  // lower right
Point(4) = { 0, -LyB, 0,lc}; // lower left

// Définition des Lignes
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Définition de la Surface

Line Loop(5) = {1,2,3,4};
Physical Line(2) = {1,3};
//Physical Line(3) = {3};
Physical Line(1) = {2,4};

Plane Surface(7) = {5};
Physical Surface(8) = {7};

Transfinite Line{1,3}=Nx; // x-resolution
Transfinite Line{2,4}=Ny; // y-resolution

Transfinite Surface(7)={1,2,3,4};

Recombine Surface(7)=0;


Mesh.Partitioner = 1; // 1 = Chaco, 2 = Metis

// Chaco
Mesh.ChacoGlobalMethod = 1;
Mesh.ChacoArchitecture = 2;

//Metis
Mesh.MetisAlgorithm = 3;
Mesh.MetisEdgeMatching = 2;
