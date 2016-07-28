/*
  Probe the flow in two locations to see if we can get a pressure trace
*/


Merge "./mesh.msh";
t = 200;
NP = NUM_PROCS;

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

// Cut the view to only the top half of the domain
Plugin(CutPlane).A = 0;
Plugin(CutPlane).B = 1;
Plugin(CutPlane).C = 0;
Plugin(CutPlane).D = 0;
Plugin(CutPlane).ExtractVolume = 1;
Plugin(CutPlane).View = v0;
Plugin(CutPlane).Run;
Delete View[v0];


// Wedge location
W = 1;
theta = ANGLE*Pi/180;  // deflection angle in radians 
xc = 5*W;              // xcoordinate of center of wedge
yc = 0;
xw = xc-0.5*W;         // position of leading point of wedge
yn = 0.5*W*Tan(theta);

// Cut plane locations btw 0.1 and 0.9 of the wedge
/*
              /|
            /  | 
          /    |
        /      |
      /        |
    /          |
  /            |
  ---x------x---

  eg. x are located at  = xc-(xc-xw)/3 and xc-(xc-xw)*2/3
  
 */

Ncuts = 10;
dx = (0.9-0.1)/(Ncuts-1);

// Loop on the cut planes
For k In {0:Ncuts-1}

// Get the cut plane coordinate
xp = xc-(xc-xw)*(0.9-k*dx);

// Do a cut plane at xp
Printf("Cut plane at x = %.6f",xp);
Plugin(CutPlane).A = 1;
Plugin(CutPlane).B = 0;
Plugin(CutPlane).C = 0;
Plugin(CutPlane).D = -xp;
Plugin(CutPlane).ExtractVolume = 0;
Plugin(CutPlane).View = v0;
Plugin(CutPlane).Run;

// Save the cut view and delete it
Save View[v0+1] Sprintf("./probe%.0f.pos",k);
Delete View[v0+1];

EndFor

Exit;
