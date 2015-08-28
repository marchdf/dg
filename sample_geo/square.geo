//
// Simple square geometry with 32^2 cells, periodic boundary conditions
// x,y = [-1,1]
//

// Points
lc = 0.3;   // length scale
Point(1) = {-1,1,0,lc};
Point(2) = {1,1,0,lc};
Point(3) = {1,-1,0,lc};
Point(4) = {-1,-1,0,lc};

// Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Surface
Line Loop(5) = {1,2,3,4};

// Periodic everywhere:
Physical Line(1) = {1,2,3,4};

// Some periodic, some reflective:
//Physical Line(3) = {1,3};
//Physical Line(1) = {2,4};

Plane Surface(7) = {5};
Physical Surface(8) = {7};

Transfinite Line{1,2,3,4}=32+1;
Transfinite Surface(7)={1,2,3,4};

Recombine Surface(7)=0;


// Multi-processor stuff
Mesh.Partitioner = 1; // 1 = Chaco, 2 = Metis

// Chaco
Mesh.ChacoGlobalMethod = 1;
Mesh.ChacoArchitecture = 2;

//Metis
Mesh.MetisAlgorithm = 3;
Mesh.MetisEdgeMatching = 2; 
