Merge "../msh1/my_line_1_dx0.msh";
Merge "./p1/dx0/rho.pos";
Merge "./p1/dx0/ux.pos";
Merge "./p1/dx0/et.pos";
Merge "./p1/dx0/p.pos";
Merge "./p1/dx0/phic.pos";
Merge "./p1/dx0/phinc.pos";
Merge "./p1/dx0/errphic.pos";
Merge "./p1/dx0/errphinc.pos";

General.Color.Background = White ;
General.Color.Foreground = Black ;
General.BackgroundGradient = 0;
General.SmallAxes = 0;
General.Trackball = 0 ;
//General.RotationX = 300;
General.RotationY = 0;
General.RotationZ = 0;
// General.ScaleX = 1.4;
// General.ScaleY = 1.5;
General.TranslationX = 0.5;
General.TranslationY = 0.8;

Geometry.Points = 0;

Mesh.Lines = 0;
Mesh.LineWidth = 1;
Mesh.Triangles = 1;
Mesh.SurfaceEdges=0;
Mesh.SurfaceFaces=0;
Mesh.ColorCarousel = 0; //Mesh coloring (0=by element type, 1=by elementary entity, 2=by physical entity, 3=by partition)

rho      = PostProcessing.NbViews-8;
ux       = rho+1;   
et       = rho+2;
p        = rho+3;
phic     = rho+4;
phinc    = rho+5;
errphic  = rho+6;
errphinc = rho+7;



// Density view
View[rho].Name = "Density";
View[rho].Axes = 2;
View[rho].Color.Axes = Black; 
View[rho].Type = 2;
View[rho].IntervalsType= 2 ; 
View[rho].RangeType = 1;
//View[rho].CustomMax = 1; 
//View[rho].CustomMin = 0; 
View[rho].LineWidth=2;
View[rho].AutoPosition = 0;
View[rho].PositionX = 85;
View[rho].PositionY = 200;
View[rho].Width = 400;
View[rho].Height = 260;
View[rho].AxesTicsX = 11;
View[rho].AxesFormatX = "%.1f";
View[rho].AxesTicsY = 10;
View[rho].SaturateValues = 1;

// Ux view
View[ux].Name = "Ux";
View[ux].Axes = 2;
View[ux].Color.Axes = Black; 
View[ux].Type = 2;
View[ux].IntervalsType= 2 ; 
View[ux].RangeType = 1;
//View[ux].CustomMax = 1; 
//View[ux].CustomMin = 0; 
View[ux].LineWidth=2;
View[ux].AutoPosition = 0;
View[ux].PositionX = 85;
View[ux].PositionY = 200;
View[ux].Width = 400;
View[ux].Height = 260;
View[ux].AxesTicsX = 11;
View[ux].AxesFormatX = "%.1f";
View[ux].AxesTicsY = 10;
View[ux].SaturateValues = 1;

// Energy view
View[et].Name = "Et";
View[et].Axes = 2;
View[et].Color.Axes = Black; 
View[et].Type = 2;
View[et].IntervalsType= 2 ; 
View[et].RangeType = 1;
//View[et].CustomMax = 1; 
//View[et].CustomMin = 0; 
View[et].LineWidth=2;
View[et].AutoPosition = 0;
View[et].PositionX = 85;
View[et].PositionY = 200;
View[et].Width = 400;
View[et].Height = 260;
View[et].AxesTicsX = 11;
View[et].AxesFormatX = "%.1f";
View[et].AxesTicsY = 10;
View[et].SaturateValues = 1;

// Pressure view
View[p].Name = "P";
View[p].Axes = 2;
View[p].Color.Axes = Black; 
View[p].Type = 2;
View[p].IntervalsType= 2 ; 
View[p].RangeType = 1;
//View[p].CustomMax = 1; 
//View[p].CustomMin = 0; 
View[p].LineWidth=2;
View[p].AutoPosition = 0;
View[p].PositionX = 85;
View[p].PositionY = 200;
View[p].Width = 400;
View[p].Height = 260;
View[p].AxesTicsX = 11;
View[p].AxesFormatX = "%.1f";
View[p].AxesTicsY = 10;
View[p].SaturateValues = 1;

// Phic view
View[phic].Name = "PhiC";
View[phic].Axes = 2;
View[phic].Color.Axes = Black; 
View[phic].Type = 2;
View[phic].IntervalsType= 2 ; 
View[phic].RangeType = 1;
//View[phic].CustomMax = 1; 
//View[phic].CustomMin = 0; 
View[phic].LineWidth=2;
View[phic].AutoPosition = 0;
View[phic].PositionX = 85;
View[phic].PositionY = 200;
View[phic].Width = 400;
View[phic].Height = 260;
View[phic].AxesTicsX = 11;
View[phic].AxesFormatX = "%.1f";
View[phic].AxesTicsY = 10;
View[phic].SaturateValues = 1;

// Phic view
View[phinc].Name = "PhiNC";
View[phinc].Axes = 2;
View[phinc].Color.Axes = Black; 
View[phinc].Type = 2;
View[phinc].IntervalsType= 2 ; 
View[phinc].RangeType = 1;
//View[phinc].CustomMax = 1; 
//View[phinc].CustomMin = 0; 
View[phinc].LineWidth=2;
View[phinc].AutoPosition = 0;
View[phinc].PositionX = 85;
View[phinc].PositionY = 200;
View[phinc].Width = 400;
View[phinc].Height = 260;
View[phinc].AxesTicsX = 11;
View[phinc].AxesFormatX = "%.1f";
View[phinc].AxesTicsY = 10;
View[phinc].SaturateValues = 1;

// Errphic view
View[errphic].Name = "ErrphiC";
View[errphic].Axes = 2;
View[errphic].Color.Axes = Black; 
View[errphic].Type = 2;
View[errphic].IntervalsType= 2 ; 
View[errphic].RangeType = 1;
//View[errphic].CustomMax = 1; 
//View[errphic].CustomMin = 0; 
View[errphic].LineWidth=2;
View[errphic].AutoPosition = 0;
View[errphic].PositionX = 85;
View[errphic].PositionY = 200;
View[errphic].Width = 400;
View[errphic].Height = 260;
View[errphic].AxesTicsX = 11;
View[errphic].AxesFormatX = "%.1f";
View[errphic].AxesTicsY = 10;
View[errphic].SaturateValues = 1;

// Errphic view
View[errphinc].Name = "ErrphiNC";
View[errphinc].Axes = 2;
View[errphinc].Color.Axes = Black; 
View[errphinc].Type = 2;
View[errphinc].IntervalsType= 2 ; 
View[errphinc].RangeType = 1;
//View[errphinc].CustomMax = 1; 
//View[errphinc].CustomMin = 0; 
View[errphinc].LineWidth=2;
View[errphinc].AutoPosition = 0;
View[errphinc].PositionX = 85;
View[errphinc].PositionY = 200;
View[errphinc].Width = 400;
View[errphinc].Height = 260;
View[errphinc].AxesTicsX = 11;
View[errphinc].AxesFormatX = "%.1f";
View[errphinc].AxesTicsY = 10;
View[errphinc].SaturateValues = 1;

// t = 0 ;

// For num In {1:400}

//   View[v0].TimeStep = t;
//   View[v1].TimeStep = t;
 
//   General.GraerrphicsWidth = 900; 
//   General.GraphicsHeight = 450;
//   General.SmallAxes = 0;
//   Sleep 0.01;		 
//   Draw;
//   Print Sprintf("bath-%03g.jpg", num);	

//   t = t+1;

// EndFor
// // // //System "mencoder 'mf://*.jpg' -mf fps=5 -o h.mpg -ovc lavc -lavcopts vcodec=mpeg1video:vhq";	
// System "mencoder mf://*.jpg -mf fps=25:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi";
