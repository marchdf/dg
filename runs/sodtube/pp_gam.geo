Merge "../msh1/my_line_1_dx4.msh";
Merge "./p1/dx4/llf/invgamma/rho.pos";
Merge "./p1/dx4/llf/invgamma/ux.pos";
Merge "./p1/dx4/llf/invgamma/et.pos";
Merge "./p1/dx4/llf/invgamma/g.pos";
Merge "./p1/dx4/llf/invgamma/p.pos";
// Merge "../msh2/my_line_2_dx4.msh";
// Merge "./p2/dx4/roe/invgamma/rho.pos";
// Merge "./p2/dx4/roe/invgamma/ux.pos";
// Merge "./p2/dx4/roe/invgamma/et.pos";
// Merge "./p2/dx4/roe/invgamma/g.pos";
// Merge "./p2/dx4/roe/invgamma/p.pos";

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

rho      = PostProcessing.NbViews-5;
ux       = rho+1;   
et       = rho+2;
g        = rho+3;
p        = rho+4;


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
s="rho.txt";
Save View[rho] s;

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

// Gamma view
View[g].Name = "G";
View[g].Axes = 2;
View[g].Color.Axes = Black; 
View[g].Type = 2;
View[g].IntervalsType= 2 ; 
View[g].RangeType = 1;
//View[g].CustomMax = 1; 
//View[g].CustomMin = 0; 
View[g].LineWidth=2;
View[g].AutoPosition = 0;
View[g].PositionX = 85;
View[g].PositionY = 200;
View[g].Width = 400;
View[g].Height = 260;
View[g].AxesTicsX = 11;
View[g].AxesFormatX = "%.1f";
View[g].AxesTicsY = 10;
View[g].SaturateValues = 1;

// Pressure view
View[p].Name = "P";
View[p].Axes = 2;
View[p].Color.Axes = Black; 
View[p].Type = 2;
View[p].IntervalsType= 2 ; 
View[p].RangeType = 2;
View[p].CustomMax = 1+1e-11; 
View[p].CustomMin = 1-1e-11; 
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


// t = 0 ;

// For num In {1:400}

//   View[v0].TimeStep = t;
//   View[v1].TimeStep = t;
 
//   General.GraphicsWidth = 900; 
//   General.GraphicsHeight = 450;
//   General.SmallAxes = 0;
//   Sleep 0.01;		 
//   Draw;
//   Print Sprintf("bath-%03g.jpg", num);	

//   t = t+1;

// EndFor
// // // //System "mencoder 'mf://*.jpg' -mf fps=5 -o h.mpg -ovc lavc -lavcopts vcodec=mpeg1video:vhq";	
// System "mencoder mf://*.jpg -mf fps=25:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi";
