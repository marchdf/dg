// Extract the growth factors based on the mass fractions Just get the
// min/max of the interface (as opposed to get_interface.geo, which
// also outputs the interface contours)

Merge "./mesh.msh";
Nt = NUM_T+1;
NP=NUM_PROCS;

// Get the isosurface
gamma = 1; 
factor = 0.50;

// Output file
fout = 'interface_minmax.dat';
Printf("# %f",Nt) > fout;

// Loop over time to write to file
For t In {0:Nt-1}

//Merge Sprintf("./y0%010.0f.pos",t);

// MPI output (or not)
If (NP == 1)
Merge Sprintf("./y0%010.0f.pos",t);
Else 
For k In {0:NP-1}
Merge Sprintf("./y0%010.0f.pos%.0f",t,k);
EndFor
EndIf

v0 = PostProcessing.NbViews-1;

Plugin(Isosurface).Value = gamma*factor ; // iso-value level
Plugin(Isosurface).ExtractVolume = 0;     // get everything at the iso-value
Plugin(Isosurface).View = v0;             // source view
Plugin(Isosurface).Run ;                  // run the plugin!
Printf("%e, %e, %e",t,View[v0+1].MinY,View[v0+1].MaxY) >> fout;

Delete View[v0+1]; // delete the isosurface view
Delete View[v0];   // delete the main view

EndFor
