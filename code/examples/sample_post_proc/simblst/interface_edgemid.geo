/*
  Extract the y position of the interface at the edge of the domain
  and at the midplane of the domain. The idea here is to track the
  interface bubble and spike. And for symmetric initial conditions of
  a sinusoidally perturbed interface, the bubble and spikes are at
  the edge and miplane of the domain (not necessarily respectively).
  
*/
   
Merge "./mesh.msh";
Nt = NUM_T+1;
NP=NUM_PROCS;

// Get the isosurface
yvalue = 0.5;
ybubble = yvalue; //0.95;
yspike  = yvalue; //0.05;

L      = 1; // length of domain
xedge  = 0;
xmid   = 0.5*L;

// Output file
fedge = 'edge.dat';
fmid  = 'mid.dat';
Printf("# %f",Nt) > fedge;
Printf("# %f",Nt) > fmid;

// Loop over time to write to file
For t In {0:Nt-1}

// MPI output (or not)
If (NP == 1)
Merge Sprintf("./y0%010.0f.pos",t);
Else 
For k In {0:NP-1}
Merge Sprintf("./y0%010.0f.pos%.0f",t,k);
EndFor
EndIf

v0 = PostProcessing.NbViews-1;


//================================================================================
// Edge vertical cut
Printf("Make a cut at the plane x = %f",xedge);
Plugin(CutPlane).A = 1 ;
Plugin(CutPlane).B = 0 ;
Plugin(CutPlane).C = 0 ;
Plugin(CutPlane).D = -xedge;
Plugin(CutPlane).View = v0;             // source view
Plugin(CutPlane).Run ;                  // run the plugin!

// extract the portion less than yvalue
Plugin(Isosurface).Value = ybubble ;       // iso-value level
Plugin(Isosurface).ExtractVolume = -1;    // get everything above the iso-value
Plugin(Isosurface).View = v0+1;           // source view
Plugin(Isosurface).Run ;                  // run the plugin!

// write to file
Printf("%e, %e",t,View[v0+2].MaxY) >> fedge;

// Clean
Delete View[v0+2]; // delete the isosurface view
Delete View[v0+1]; // delete the cut view

//================================================================================
// Midplane vertical cut
Printf("Make a cut at the plane x = %f",xmid);
Plugin(CutPlane).A = 1 ;
Plugin(CutPlane).B = 0 ;
Plugin(CutPlane).C = 0 ;
Plugin(CutPlane).D = -xmid;
Plugin(CutPlane).View = v0;             // source view
Plugin(CutPlane).Run ;                  // run the plugin!

// extract the portion less than yvalue
Plugin(Isosurface).Value = yspike ;       // iso-value level
Plugin(Isosurface).ExtractVolume = 1;     // get everything above the iso-value
Plugin(Isosurface).View = v0+1;           // source view
Plugin(Isosurface).Run ;                  // run the plugin!

// write to file
Printf("%e, %e",t,View[v0+2].MinY) >> fmid;

// Clean
Delete View[v0+2]; // delete the isosurface view
Delete View[v0+1]; // delete the cut view

Delete View[v0];
EndFor
