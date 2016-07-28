/*
  Output the mach contours for the bubbly flow over a wedge problem
*/


Merge "./mesh.msh";
t = 100;
NP = 1;

General.Color.Background = Black ;
General.Color.Foreground = White ;
General.Color.Text = White ;
General.BackgroundGradient = 0;
General.SmallAxes = 0;
General.Trackball = 0 ;
General.RotationX = 0;
General.RotationY = 0;
General.RotationZ = 00;
General.ScaleX = 9;
General.ScaleY = 9;
General.TranslationX = -5;
General.TranslationY = 0;
General.GraphicsWidth = 2400 ;
General.GraphicsHeight = 1280 ;
General.GraphicsFontSize = 24;
General.GraphicsFontSizeTitle = 24;

Geometry.Points = 0;

Mesh.Lines = 0;
Mesh.LineWidth = 2;
Mesh.Triangles = 1;
Mesh.SurfaceEdges=0;
Mesh.SurfaceFaces=0;
Mesh.ColorCarousel = 0; //Mesh coloring (0=by element type, 1=by elementary entity, 2=by physical entity, 3=by partition)

PostProcessing.HorizontalScales=0;

// Load the final time step
If (NP == 1)
Merge Sprintf("./rho%010.0f.pos",t);
Merge Sprintf("./ux%010.0f.pos",t);   // View 1
Merge Sprintf("./uy%010.0f.pos",t);   // View 2
Merge Sprintf("./p%010.0f.pos",t);    // View 3
Merge Sprintf("./pinf%010.0f.pos",t); // View 4
Merge Sprintf("./g%010.0f.pos",t);    // View 5
EndIf
If (NP != 1)
For k In {0:NP-1}
Merge Sprintf("./rho%010.0f.pos%.0f",t,k);  // View 0
EndFor 
For k In {0:NP-1}
Merge Sprintf("./ux%010.0f.pos%.0f",t,k);   // View 1
EndFor 
For k In {0:NP-1}
Merge Sprintf("./uy%010.0f.pos%.0f",t,k);   // View 2
EndFor 
For k In {0:NP-1}
Merge Sprintf("./p%010.0f.pos%.0f",t,k);    // View 3
EndFor 
For k In {0:NP-1}
Merge Sprintf("./pinf%010.0f.pos%.0f",t,k); // View 4
EndFor
For k In {0:NP-1}
Merge Sprintf("./g%010.0f.pos%.0f",t,k);    // View 5
EndFor
EndIf

v0 = 0;

// Extract the time step for each one
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).TimeStep = t;
Plugin(MathEval).View = 0;
Plugin(MathEval).Run;

// Extract the time step for each one
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).TimeStep = t;
Plugin(MathEval).View = 1;
Plugin(MathEval).Run;

// Extract the time step for each one
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).TimeStep = t;
Plugin(MathEval).View = 2;
Plugin(MathEval).Run;

// Extract the time step for each one
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).TimeStep = t;
Plugin(MathEval).View = 3;
Plugin(MathEval).Run;

// Extract the time step for each one
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).TimeStep = t;
Plugin(MathEval).View = 4;
Plugin(MathEval).Run;

// Extract the time step for each one
Plugin(MathEval).Expression0 = "v0";
Plugin(MathEval).TimeStep = t;
Plugin(MathEval).View = 5;
Plugin(MathEval).Run;

// remove the old views
Delete View[v0];
Delete View[v0];
Delete View[v0];
Delete View[v0];
Delete View[v0];
Delete View[v0];

// Take the magnitude of the velocity
Plugin(MathEval).Expression0 = "Sqrt(v0*v0+w0*w0)";
Plugin(MathEval).TimeStep = 0;
Plugin(MathEval).View = 1;
Plugin(MathEval).OtherTimeStep = 0;
Plugin(MathEval).OtherView = 2;
Plugin(MathEval).Run;

Delete View[2];
Delete View[1];
// rho -> view 0
// p   -> view 1
// pinf-> view 2
// g   -> view 3
// |u| -> view 4


// Calculate the pressure sum (p+pinf)
Plugin(MathEval).Expression0 = "v0+w0";
Plugin(MathEval).TimeStep = 0;
Plugin(MathEval).View = 1;
Plugin(MathEval).OtherTimeStep = 0;
Plugin(MathEval).OtherView = 2;
Plugin(MathEval).Run;

Delete View[2];
Delete View[1];
// rho    -> view 0
// g      -> view 1
// |u|    -> view 2
// p+pinf -> view 3

// Calculate g*(p+pinf)
Plugin(MathEval).Expression0 = "v0*w0";
Plugin(MathEval).TimeStep = 0;
Plugin(MathEval).View = 1;
Plugin(MathEval).OtherTimeStep = 0;
Plugin(MathEval).OtherView = 3;
Plugin(MathEval).Run;

Delete View[3];
Delete View[1];
// rho        -> view 0
// |u|        -> view 1
// g*(p+pinf) -> view 2

// Calculate the local speed of sound
Plugin(MathEval).Expression0 = "Sqrt(v0/w0)";
Plugin(MathEval).TimeStep = 0;
Plugin(MathEval).View = 2;
Plugin(MathEval).OtherTimeStep = 0;
Plugin(MathEval).OtherView = 0;
Plugin(MathEval).Run;

Delete View[2];
Delete View[0];
// |u|                  -> view 0
// sqrt(g*(p+pinf)/rho) -> view 1


// Calculate the Mach number
Plugin(MathEval).Expression0 = "v0/w0";
Plugin(MathEval).TimeStep = 0;
Plugin(MathEval).View = 0;
Plugin(MathEval).OtherTimeStep = 0;
Plugin(MathEval).OtherView = 1;
Plugin(MathEval).Run;

Delete View[1];
Delete View[0];
// u/sqrt(g*(p+pinf)/rho) -> view 0

// Custom Mach number
View[v0].Name = "Mach";
View[v0].Light = 0;
//View[v0].ScaleType = 2; 
View[v0].RangeType = 2;
View[v0].CustomMax = 3;
View[v0].CustomMin = 0;
View[v0].SaturateValues = 1;
//View[v0].Format = "%5.5g";
View[v0].NbIso = 4;
View[v0].ShowTime = 0;
View[v0].Color.Axes = White;
View[v0].OffsetX = 0.0;
View[v0].OffsetY = 0.0;
View[v0].AutoPosition = 0;
View[v0].PositionX = 1270;
View[v0].PositionY = 1000;
View[v0].Width = 20;
View[v0].Height = 250;
View[v0].Axes = 1;
View[v0].AxesFormatY = "%.2g";
View[v0].AxesTicsX = 2;
View[v0].AxesTicsY = 5;
View[v0].ColormapNumber = 5;
View[v0].ColormapInvert = 1;
View[v0].ColormapSwap = 1;
View[v0].Clip = 2;

Draw;
Print.Background = 1;
Print Sprintf("./mach.jpg");

// Exit;


// // // Custom the density view
// // View[v0].Light=0;
// // View[v0].Visible = 0;
// // // Adapt visualization
// // View[v0].AdaptVisualizationGrid = 1;
// // View[v0].MaxRecursionLevel = 1;
// // View[v0].TargetError = -1; 

// // // Take the gradient of density
// // Plugin(Gradient).View = v0;
// // Plugin(Gradient).Run ;                  // run the plugin!
// // View[v0+1].Visible = 0;

// // // Take the norm of the gradient of density
// // Plugin(MathEval).Expression0 = "Sqrt(v0^2+v1^2+v2^2)";
// // Plugin(MathEval).TimeStep = 0;
// // Plugin(MathEval).View = v0+1;
// // Plugin(MathEval).Run;

// // // Find the Maximum of the norm of the gradient
// // Plugin(MinMax).View = v0+2;
// // Plugin(MinMax).Run;
// // Max = View[v0+4].Max ;

// // // Take the norm of the gradient of density
// // Plugin(MathEval).Expression0 = Sprintf("Exp(-20*Sqrt(v0^2+v1^2+v2^2)/%e)",Max);
// // Plugin(MathEval).TimeStep = 0;
// // Plugin(MathEval).View = v0+1;
// // Plugin(MathEval).Run;

// // Delete View[v0+2]; // delete the norm view
// // Delete View[v0+2]; // delete the min view
// // Delete View[v0+2]; // delete the max view

// // // Custom the Schlieren view
// // View[v0+2].Light=0;
// // View[v0+2].Name = "Sch.";
// // View[v0+2].ShowScale = 0;
// // View[v0+2].RangeType = 2;
// // View[v0+2].CustomMax = 1;
// // View[v0+2].CustomMin = 0.4;
// // View[v0+2].OffsetX = 0.04;
// // View[v0+2].SaturateValues = 1;
// // View[v0+2].Format = "%5.5g";
// // View[v0+2].ShowTime = 2;
// // View[v0+2].Color.Axes = White;
// // View[v0+2].ColormapNumber = 9;
// // View[v0+2].AutoPosition = 0;
// // View[v0+2].PositionX = 850;
// // View[v0+2].PositionY = 150;
// // View[v0+2].Width = 20;
// // View[v0+2].Height = 430;
// // View[v0+2].Axes = 1;
// // View[v0+2].AxesFormatY = "%.2g";
// // View[v0+2].AxesTicsX = 0;
// // View[v0+2].AxesTicsY = 3;

// // // // Adapt visualization
// // // View[v0+2].AdaptVisualizationGrid = 1;
// // // View[v0+2].MaxRecursionLevel = 1;
// // // View[v0+2].TargetError = -1; 


