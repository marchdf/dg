/*
  Output a picture of the wedge, later used to measure the angles.  
*/


Merge "./mesh.msh";
t = 200;
NP = NUM_PROCS;

General.Color.Background = Black ;
General.Color.Foreground = White ;
General.Color.Text = White ;
General.BackgroundGradient = 0;
General.SmallAxes = 0;
General.Trackball = 0 ;
General.RotationX = 0;
General.RotationY = 0;
General.RotationZ = 0;
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

// view options
View[v0].Light=0;
View[v0].ColormapNumber = 5;
View[v0].ColormapInvert = 1;
View[v0].ColormapSwap = 1;
View[v0].ShowScale = 0;
View[v0].Light = 0;

// Adapt visualization
View[v0].AdaptVisualizationGrid = 1;
View[v0].MaxRecursionLevel = 1;
View[v0].TargetError = -1; 

Draw;
Print.Background = 1;
Print Sprintf("./picture.jpg");

// print again but with a different zoom
PostProcessing.HorizontalScales=0;
View[v0].NbIso = 4;
View[v0].ShowTime = 0;
View[v0].ShowScale = 1;
View[v0].AutoPosition = 0;
View[v0].PositionX = 110;
View[v0].PositionY = 250;
View[v0].Width = 20;
View[v0].Height = 500;
General.ScaleX = 3.5;
General.ScaleY = 3.5;
Draw;
Print.Background = 1;
Print Sprintf("./picture_unzoom.jpg");

Exit;

