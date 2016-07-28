//
// Load pressure density data (used by viz.geo).
// This one was used for the drop hitting a wall problem.
// 

Merge Sprintf("./rho%010.0f.pos",num);
Merge Sprintf("./p%010.0f.pos",num);
// For k In {0:NP-1}
// Merge Sprintf("./rho%010.0f.pos%.0f",num,k);
// EndFor 


// Assign values to the views
v0 = PostProcessing.NbViews-2;
v1 = PostProcessing.NbViews-1;

//
// Custom density
//

// Set the clipping planes for density
General.Clip0A = 0;
General.Clip0B = 1;
General.Clip0C = 0;
General.Clip0D = 0; 

General.Clip1A = 1;
General.Clip1B = 0;
General.Clip1C = 0;
General.Clip1D = 8; 

General.Clip2A = 0;
General.Clip2B = -1;
General.Clip2C = 0;
General.Clip2D = 6; 

// Options
View[v0].Name = "Density";
View[v0].ScaleType = 2; 
View[v0].RangeType = 2;
View[v0].CustomMax = 1000;
View[v0].CustomMin = 0.1;
View[v0].SaturateValues = 1;
//View[v0].Format = "%5.5g";
View[v0].NbIso = 4;
View[v0].ShowTime = 0;
View[v0].Color.Axes = White;
View[v0].OffsetX = 0.0;
View[v0].OffsetY = 0.0;
View[v0].AutoPosition = 0;
View[v0].PositionX = 1080;
View[v0].PositionY = -220;
View[v0].Width = 20;
View[v0].Height = 330;
View[v0].Axes = 1;
View[v0].AxesFormatY = "%.2g";
View[v0].AxesTicsX = 2;
View[v0].AxesTicsY = 5;
View[v0].ColormapNumber = 5;
View[v0].ColormapInvert = 1;
View[v0].ColormapSwap = 1;
View[v0].Clip = 1+2+4;
View[v0].Light = 0;

//Include "map.map";

//
// Custom pressure
//

// Set the clipping planes for pressure
General.Clip3A = 0;
General.Clip3B = -1;
General.Clip3C = 0;
General.Clip3D = 0; 

General.Clip4A = 0;
General.Clip4B = 1;
General.Clip4C = 0;
General.Clip4D = 6; 

// Options
View[v1].Name = "Pressure";
//View[v1].ScaleType = 2; 
View[v1].RangeType = 3;//2;
View[v1].CustomMax = 20;
View[v1].CustomMin = -5;
View[v1].SaturateValues = 1;
//View[v1].Format = "%5.5g";
View[v1].NbIso = 4;
View[v1].ShowTime = 0;
View[v1].Color.Axes = White;
View[v1].OffsetX = 0.0;
View[v1].OffsetY = 0.0;
View[v1].AutoPosition = 0;
View[v1].PositionX = 1080;
View[v1].PositionY = 230;
View[v1].Width = 20;
View[v1].Height = 330;
View[v1].Axes = 1;
View[v1].AxesFormatY = "%.2g";
View[v1].AxesTicsX = 2;
View[v1].AxesTicsY = 5;
View[v1].ColormapNumber = 5;
View[v1].ColormapInvert = 1;
View[v1].ColormapSwap = 1;
View[v1].Clip = 2+8+16;
View[v1].Light = 0;
//View[v1].TransformYY = -1;




// // Adapt visualization
// View[v0].AdaptVisualizationGrid = 1;
// View[v0].MaxRecursionLevel = 2;
// View[v0].TargetError = -1;


T = t*Dtout;

// Add a title
Plugin(Annotate).Text = Sprintf("t = %5.2f",T);
Plugin(Annotate).X = 1.e5; // By convention, for window coordinates a value greater than 99999 represents the center.
Plugin(Annotate).Y = 70 ;
Plugin(Annotate).FontSize = 24;
Plugin(Annotate).Align = "Center" ;
Plugin(Annotate).View = 0;
Plugin(Annotate).Run ;
View[2].Color.Text2D = {0,0,0}; //{255,255,255};

//
// Overlay mass fraction contours (derived from pinf)
//
// View 2: load pinf
Merge Sprintf("./pinf%010.0f.pos",num);
// For k In {0:NP-1}
// Merge Sprintf("./pinf%010.0f.pos%.0f",num,k);
// EndFor
Max = View[3].Max;

// View 3: normalized pinf
Plugin(MathEval).Expression0 = Sprintf("v0/(%e)",Max);
Plugin(MathEval).TimeStep = -1;
Plugin(MathEval).View = 3;
Plugin(MathEval).Run;
Delete View[3]; // view 3 is now view 2

//View[3].Name = "Pressure";
View[3].IntervalsType = 1;
View[3].NbIso = 100;
View[3].RangeType = 2;
View[3].CustomMax = 0.51;
View[3].CustomMin = 0.49;
View[3].SaturateValues = 1;
View[3].ShowScale = 0;
View[3].LineWidth = 4;
View[3].Light = 0;
View[3].ColorTable = {{0,0,0,255}};//{{255, 255, 255,255}};//{{238, 46, 47,255}};
View[3].Clip = 2+8+16;
View[3].Light = 0;

Sleep 0.0001;
Draw;
Print.Background = 1;
Print Sprintf("./pics_rho/sd%03g.jpg", num);

Delete View[3];
Delete View[2];
Delete View[1];
Delete View[0];
