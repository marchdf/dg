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
scale = w/NUM_CELLS;
Point(1)={0,0,0,scale};
Point(2)={w,0,0,scale};
Point(3)={w,h,0,scale};
Point(4)={0,h,0,scale};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {3,4,1,2};
Physical Line(2) = {1,2,3,4};
Plane Surface(6) = {5};
Physical Surface(1) = {6};
