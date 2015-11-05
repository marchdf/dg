/* Mesh for the vortex transport problem defined in the High-Order workshop (2014)
   (1)------(2)
    |        |
    |        |
    |        |
   (4)------(3)
*/

Lx = 1;
Ly = 1;
N = RESOLUTION;  // number of cells

// Points
lc = 0.3;   // length scale
Point(1) = {-Lx, Ly,0,lc};
Point(2) = { Lx, Ly,0,lc};
Point(3) = { Lx,-Ly,0,lc};
Point(4) = {-Lx,-Ly,0,lc};

// Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Surface
Line Loop(5) = {1,2,3,4};

// Periodic everywhere:
Physical Line(1) = {1,2,3,4};

Plane Surface(7) = {5};
Physical Surface(8) = {7};

Transfinite Line{1,2,3,4}=N+1;
Transfinite Surface(7)={1,2,3,4};

//Recombine Surface(7)=0;

Mesh.Partitioner = 1; // 1 = Chaco, 2 = Metis

// Chaco
Mesh.ChacoGlobalMethod = 1;
Mesh.ChacoArchitecture = 2;

//Metis
Mesh.MetisAlgorithm = 3;
Mesh.MetisEdgeMatching = 2; 
