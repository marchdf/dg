// Just a 1D sample mesh with many elements

lc = 0.3;

// Domain length
Lx = 1;
//LyA = 20*Lx;
//LyB = 30*Lx;
LyB = 2;
LyA = 3;

Nx = 100+1;

// Point definitions
Point(1) = {-LyB,0,0,lc};
Point(2) = { LyA,0,0,lc};

// Line definitions
Line(1) = {1,2};

// Periodic everywhere:
Physical Point(1) = {1,2}; //periodic BC

// Surface definitions
//Physical Point(2) = {1,2}; //tells Marc's code to enforce zero-gradient BC
Physical Line(1) = {1,2};

Transfinite Line{1}=Nx;
