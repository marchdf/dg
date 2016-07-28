/* Output the circulation and the density weighted enstrophy to a file.

   This is the same as circ_enst_half.geo, EXCEPT it gets the
   circulation in the whole domain.

 */

Merge "./mesh.msh";
Nt = NUM_T+1;
NP=NUM_PROCS;

// Loop on time and output the values
fout = 'circ_enst.dat';
Printf("# Circulation and enstrophy") > fout;

For t In {0:Nt-1}

// MPI output
If (NP == 1)
Merge Sprintf("./rho%010.0f.pos",t);
Merge Sprintf("./ux%010.0f.pos",t);
Merge Sprintf("./uy%010.0f.pos",t);
Else 
For k In {0:NP-1}
Merge Sprintf("./rho%010.0f.pos%.0f",t,k);
Merge Sprintf("./ux%010.0f.pos%.0f",t,k);
Merge Sprintf("./uy%010.0f.pos%.0f",t,k);
EndFor
EndIf

v0 = PostProcessing.NbViews-3;
v1 = PostProcessing.NbViews-2;
v2 = PostProcessing.NbViews-1;

// Extract the right times, delete the old ones
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).View = v0;
Plugin(MathEval).TimeStep = t;
Plugin(MathEval).Run;
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).View = v1;
Plugin(MathEval).TimeStep = t;
Plugin(MathEval).Run;
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).View = v2;
Plugin(MathEval).TimeStep = t;
Plugin(MathEval).Run;
Delete View[2];
Delete View[1];
Delete View[0];

// View 3: vector view for velocity
Plugin(Scal2Vec).NameNewView="U";
Plugin(Scal2Vec).ViewX=v1;
Plugin(Scal2Vec).ViewY=v2;
Plugin(Scal2Vec).Run;

// View 4: Take the curl of velocity = vorticity
Plugin(Curl).View = 3;
Plugin(Curl).Run;

// View 5: Integrate vorticity to get circulation
Plugin(Integrate).View = 4;
Plugin(Integrate).Run;

// // View 6: Get the enstrophy (multiply the half domain one by 2 to get it in the whole domain)
Plugin(MathEval).Expression0 = "2*v0*(w0*w0+w1*w1+w2*w2)";
Plugin(MathEval).View = v0;
Plugin(MathEval).OtherView = 4;
Plugin(MathEval).TimeStep = 0;
Plugin(MathEval).Run;

// View 7: Integrate previous to get enstrophy
Plugin(Integrate).View = 6;
Plugin(Integrate).Run;

// View 8: Extract number circulation view
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).View = 5;
Plugin(MathEval).OtherView = -1;
Plugin(MathEval).Run;

// View 9: Extract number enstrophy view
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).View = 7;
Plugin(MathEval).OtherView = -1;
Plugin(MathEval).Run;

Printf("%f, %e, %e",t,View[8].Max,View[9].Max) >> fout;

Delete View[9];
Delete View[8];
Delete View[7];
Delete View[6];
Delete View[5];
Delete View[4];
Delete View[3];
Delete View[2];
Delete View[1];
Delete View[0];

EndFor
