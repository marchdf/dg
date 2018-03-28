//This file taken shamelessly from http://matveichev.blogspot.com/2013/12/building-hexagonal-meshes-with-gmsh.html
//and modified a bit for my purposes.

// Syntax.: Point(n) = {x, y, z, hs};
// As we'll use transfinite algorithm hs can be arbitrary
// or even omitted from mesh description

Ld = 1.0*Pi;
Point(1) = {-Ld, -Ld, -Ld, 1.0}; //--- 
Point(2) = { Ld, -Ld, -Ld, 1.0}; //+--
Point(3) = {-Ld,  Ld, -Ld, 1.0}; //-+-
Point(4) = { Ld,  Ld, -Ld, 1.0}; //++-
Point(5) = { Ld,  Ld,  Ld, 1.0}; //+++
Point(6) = { Ld, -Ld,  Ld, 1.0}; //+-+
Point(7) = {-Ld,  Ld,  Ld, 1.0}; //-++
Point(8) = {-Ld, -Ld,  Ld, 1.0}; //--+

//connect them with lines

Line(1) = {3, 7}; //-+- to -++
Line(2) = {7, 5}; //-++ to +++
Line(3) = {5, 4}; //+++ to ++-
Line(4) = {4, 3}; //++- to -+-
Line(5) = {3, 1}; //-+- to ---
Line(6) = {2, 4}; //+-- to ++-
Line(7) = {2, 6}; //+-- to +-+
Line(8) = {6, 8}; //+-+ to --+
Line(9) = {8, 1}; //--+ to ---
Line(10) = {1, 2}; //--- to +--
Line(11) = {8, 7}; //--+ to -++
Line(12) = {6, 5}; //+-+ to +++

//Physical Line(80) = {1,2,3,4,5,6,7,8,9,10,11,12};

Line Loop(13) = {7, 8, 9, 10};
Plane Surface(14) = {13};
Line Loop(15) = {6, 4, 5, 10};
Plane Surface(16) = {15};
Line Loop(17) = {3, 4, 1, 2};
Plane Surface(18) = {17};
Line Loop(19) = {12, -2, -11, -8};
Plane Surface(20) = {19};
Line Loop(21) = {7, 12, 3, -6};
Plane Surface(22) = {21};
Line Loop(23) = {9, -5, 1, -11};
Plane Surface(24) = {23};

// Periodic everywhere:
//
//Physical Surface(80) = {13,14};

Surface Loop(25) = {14, 22, 20, 18, 16, 24};

//Physical Surface(80) = {14,16,18,20,22,24};

//Physical Line(80) = {6};

//and volume

Volume(26) = {25};



//And the next lines will convert default generation strategy from tetrahedral meshes to hexagonal

OverRes = 21;
Nx = OverRes+1; //resolution in x direction
Ny = OverRes+1; //resolution in y direction
Nz = OverRes+1; //resolution in z direction
//Transfinite Line{1,2,3,4} = Nx; //this is where mesh resolution is defined
//Transfinite Surface(7)={1,2,3,4};
Transfinite Line{2,4,8,10} = Nx;
Transfinite Line{5,6,11,12} = Ny;
Transfinite Line{1,3,7,9} = Nz;
//Transfinite Line "*" = 100 Using Bump 0.25;
Transfinite Surface "*";

Recombine Surface "*";
//Physical Surface(80) = {14,16,18,20,22,24};
Transfinite Volume "*";
//Add physical tag to surfaces for periodic BC:
Physical Surface(1) = {14,16,18,20,22,24};
//Need a physical tag for the mesh volume so the elements make it in to the .msh
Physical Volume(8) = {26};

//As the lengths of the cube sides are equal, densities of the 1D mesh on every side are equal (also I've used keyword Bump to grade the mesh near the corners). To use generated mesh in OpenFOAM one has to define physical groups also:

//Physical Surface("top") = {20};
//Physical Surface("bottom") = {16};
//Physical Surface("sides") = {14, 24, 18, 22};
//Physical Volume("cube") = {26};

//Multi-processor stuff
Mesh.Partitioner = 1; // 1 = Chaco, 2 = Metis

// Chaco
Mesh.ChacoGlobalMethod = 1;
Mesh.ChacoArchitecture = 2;

//Metis
Mesh.MetisAlgorithm = 3;
Mesh.MetisEdgeMatching = 2; 


