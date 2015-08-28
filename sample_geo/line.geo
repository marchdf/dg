//
// 1D line with zero-gradient boundary conditions
// x = [0,10] and 32 cells
//


// Number of points
NUM=32+1; // i.e. 32 cells

// characteristic length scale
lc = 0.3;   

// Point definition
Point(1) = {0,0,0,lc};
Point(2) = {10,0,0,lc};

// Line definition
Line(1) = {1,2};

// Surface definition
Physical Point(2) = {1,2};
Physical Line(1) = {1,2,3,4};

// Regular mesh
Transfinite Line{1}=NUM;

