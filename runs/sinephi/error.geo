
phi        = 0;
phi_diff2  = phi+1;
phi_error  = phi+2;

For pdir In {0:4}

For dxdir In {0:4}

Merge Sprintf("../msh%01g/my_line_%01g_dx%01g.msh", pdir, pdir, dxdir);
Merge Sprintf("./p%01g/dx%01g/phic.pos", pdir, dxdir);

Plugin(MathEval).Expression0 = "(w0-v0)*(w0-v0)";
Plugin(MathEval).View = phi;
Plugin(MathEval).TimeStep = 0;
Plugin(MathEval).OtherView = phi;
Plugin(MathEval).OtherTimeStep = 200;
Plugin(MathEval).Run;
Plugin(Integrate).View=phi_diff2;
Plugin(Integrate).Run;

Save View[phi_error] Sprintf("./p%01g/dx%01g/phic_error.dat", pdir, dxdir);

Delete View[phi_error];
Delete View[phi_diff2];
Delete View[phi];

////////////////////////////////////////////////////////////
Merge Sprintf("./p%01g/dx%01g/phinc.pos", pdir, dxdir);

Plugin(MathEval).Expression0 = "(w0-v0)*(w0-v0)";
Plugin(MathEval).View = phi;
Plugin(MathEval).TimeStep = 0;
Plugin(MathEval).OtherView = phi;
Plugin(MathEval).OtherTimeStep = 200;
Plugin(MathEval).Run;
Plugin(Integrate).View=phi_diff2;
Plugin(Integrate).Run;

Save View[phi_error] Sprintf("./p%01g/dx%01g/phinc_error.dat", pdir, dxdir);

Delete View[phi_error];
Delete View[phi_diff2];
Delete View[phi];

EndFor

EndFor
