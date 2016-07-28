// Extract the growth factors based on the mass fractions. It does two
// things: a) save the min/max of the interface to a file and b) save
// the isosurface of the interface as a .pos file

Merge "./mesh.msh";
Nt = NUM_T+1;
NP=NUM_PROCS;

// Get the isosurface
yvalue = 0.5;

// Output file for the min and max
fout = 'interface_minmax.dat';
Printf("# %f",Nt) > fout;

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

Plugin(Isosurface).Value = yvalue ;       // iso-value level
Plugin(Isosurface).ExtractVolume = 0;     // get everything at the iso-value
Plugin(Isosurface).View = v0;             // source view
Plugin(Isosurface).Run ;                  // run the plugin!

// Save the min/max of the interface
Printf("%e, %e, %e",t,View[v0+1].MinY,View[v0+1].MaxY) >> fout;

// Save the interface itself
Save View[v0+1] Sprintf("./interface%03.0f.pos",t); 

Delete View[v0+1]; // delete the isosurface view
Delete View[v0];   // delete the main view

EndFor
