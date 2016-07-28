/* Calculate the pressure on the wedge (behind the shock).
   You have to provide (e.g. via sed) the angle of the wedge.
*/


Merge "./mesh.msh";
t = 200;
NP = NUM_PROCS;

// Loop on time and output the values
fout = 'pwedge.dat';
Printf("# Pressure on the wedge: xp, yp, pressure") > fout;

// Load the final time step
If (NP == 1)
Merge Sprintf("./p%010.0f.pos",t);
EndIf
If (NP != 1)
For k In {0:NP-1}
Merge Sprintf("./p%010.0f.pos%.0f",t,k);
EndFor 
EndIf
v0 = PostProcessing.NbViews-1;

// Extract the final time step
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).TimeStep = t;
Plugin(MathEval).View = v0;
Plugin(MathEval).Run;
Delete View[v0];

// Wedge location
W = 1;
theta = ANGLE*Pi/180;  // deflection angle in radians 
xc = 5*W;              // xcoordinate of center of wedge
yc = 0;
xw = xc-0.5*W;         // position of leading point of wedge
yn = 0.5*W*Tan(theta);

// probe location (halfway up the wedge)
/*
              /|
            /  | 
          /    |
        x      |
      /        |
    /          |
  /            |
  --------------

  x is located at (xp,yp) = (xc-0.5*(xc-xw),yc+0.5*(yn-yc))
  
 */

xp = xc-0.5*(xc-xw);
yp = yc+0.5*(yn-yc);

// Probe the pressure on the wedge, behind the shock
Printf("Probe at (xp,yp)= (%.6f,%.6f)",xp,yp);
Plugin(Probe).X = xp;
Plugin(Probe).Y = yp;
Plugin(Probe).View = v0;
Plugin(Probe).Run;

Printf("%.16f, %.16f, %.16f",xp,yp,View[v0+1].Max) >> fout;
