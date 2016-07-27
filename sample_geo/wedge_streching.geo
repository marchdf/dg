// Mesh for two phase flow over a wedge
// Generates a mesh for the whole domain
/*

  (1)-------(2)--------------------(3)--------(4)
   |         |                      |          |
   |         |                      |          |
   |         |                      |          |
   |         |         (13)         |          |
   |         |        /    \        |          |
   |         |       /      \       |          |
  (10)------(11)---(16)----(14)----(12)-------(5)
   |         |       \      /       |          |
   |         |        \    /        |          |
   |         |         (15)         |          |
   |         |                      |          |
   |         |                      |          |
   |         |                      |          |
  (9)-------(8)--------------------(7)--------(6)

  in the paper (Eddington 1967):
  - the wedge is 0.5 in wide
  - the test section height is 1.6 in
 */

// dimentions
w = 0.5; // inches
h = 1.6; // inches
L_ND = w;

// Dimensionless dimensions
W = w/L_ND;            // width of wedge
H = 0.5*h/L_ND;        // half height of channel
theta = 10.0*Pi/180.0; // deflection angle in radians
xc = 5*W;              // xcoordinate of center of wedge
char_mesh_size = 0.01*W;
char_mesh_size_large = 10*char_mesh_size;

// for 4 degree angle:
Mesh.RandomFactor = 1e-8;

// Derived quantities
Lx = 5*H;       // length in x of higher resolution portion
Lxl = 2*H;      // length in x of grid stretching portion
xw = xc-0.5*W;  // position of leading point of wedge
xe = xc+0.5*W;  // position of trailing point of wedge
yn = 0.5*W*Tan(theta);
ys = -yn;


// Définition des Points
lc = char_mesh_size;   // longueur caractéristique
lc_large = char_mesh_size_large;   // longueur caractéristique
Point(1) = {-Lxl, H,0,lc_large};
Point(2) = { 0, H,0,lc};  
Point(3) = { Lx, H,0,lc};
Point(4) = { Lx+Lxl,H,0,lc_large};
Point(5) = { Lx+Lxl,0,0,lc_large}; 
Point(6) = { Lx+Lxl,-H,0,lc_large};
Point(7) = { Lx, -H,0,lc};
Point(8) = { 0, -H,0,lc};
Point(9) = {-Lxl, -H,0,lc_large};
Point(10)= {-Lxl,0,0,lc_large};
Point(11)= {0,0,0,lc};
Point(12)= {Lx,0,0,lc};
Point(13)={xc,yn,0,lc};
Point(14) = {xe,0,0,lc};
Point(15) = {xc,ys,0,lc};
Point(16) = {xw,0,0,lc};

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

Line(11) = {3,12};
Line(12) = {12,7};

Line(13) = {8,11};
Line(14) = {11,2};

Line(15) = {13,14};
Line(16) = {14,15};
Line(17) = {15,16};
Line(18) = {16,13};

// Définition de la Surface
Line Loop(44) = {1,-14,-13,8,9,10}; // grid stretch region left
Line Loop(45) = {2,11,12,7,13,14};  // high resolution region
Line Loop(46) = {3,4,5,6,-12,-11};  // grid stretch region right
Line Loop(47) = {15,16,17,18};      // wedge
Physical Line(2) = {4,5,9,10};
Physical Line(3) = {1,2,3,6,7,8,15,16,17,18};


Plane Surface(80) = {45,47}; // this creates a hole in the mesh
Plane Surface(81) = {44};
Plane Surface(82) = {46}; 
Physical Surface(90) = {80,81,82};



Mesh.Partitioner = 1; // 1 = Chaco, 2 = Metis

// Chaco
Mesh.ChacoGlobalMethod = 1;
Mesh.ChacoArchitecture = 2;

//Metis
Mesh.MetisAlgorithm = 3;
Mesh.MetisEdgeMatching = 2; 
