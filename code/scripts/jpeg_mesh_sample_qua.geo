// Same as jpeg_mesh_sample_tri.geo but regular cartesian grid
//
// This meshes a jpg image. And saves the image data in the GMSH mesh
// format so that it can be read by the code.
// Modified by original was downloaded from http://geuz.org/photos/cg/cg.geo

// Merge the image (this will create a new post-processing view,
// View[0])
Merge "IMAGENAME";

// Modify the normalized pixel values to make them reasonnable
// characteristic lengths for the mesh
Plugin(ModifyComponent).Expression = "v0 * 100";
Plugin(ModifyComponent).Run;

// Apply the view as the current background mesh This will mesh with
// the picture (refine where refinement is needed, ie lots of picture
// details)). Leave commented if you want a uniform mesh
//Background Mesh View[0];

// Build a simple geometry on top of the background mesh (or jpeg)
w = View[0].MaxX;
h = View[0].MaxY;
scale = w;
Point(1)={0,0,0,scale};
Point(2)={w,0,0,scale};
Point(3)={w,h,0,scale};
Point(4)={0,h,0,scale};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Surface
Line Loop(5) = {1,2,3,4};

// zero gradient everywhere
Physical Line(2) = {1,2,3,4};

Plane Surface(6) = {5};
Physical Surface(1) = {6};

// regular grid
Transfinite Line{1,3}=NUM_CELLS+1;
Transfinite Line{2,4}=NUM_CELLS*h/w+1;
Transfinite Surface(6)={1,2,3,4};
Recombine Surface(6)=0;
