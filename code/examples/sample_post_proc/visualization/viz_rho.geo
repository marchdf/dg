//
// View density data (to be used by viz.geo)
// 

Merge Sprintf("./rho%010.0f.pos",num);
// For k In {0:NP-1}
// Merge Sprintf("./rho%010.0f.pos%.0f",num,k);
// EndFor 

// Assign values to the views
v0 = PostProcessing.NbViews-1;

// Custom density
View[v0].Name = "Density";
View[v0].ScaleType = 1; 
View[v0].RangeType = 2;
View[v0].CustomMax = 3.1;
View[v0].CustomMin = 0.9;
View[v0].SaturateValues = 1;
//View[v0].Format = "%5.5g";
View[v0].NbIso = 4;
View[v0].ShowTime = 0;
View[v0].Color.Axes = White;
View[v0].OffsetX = 0.0;
View[v0].OffsetY = 0.0;
View[v0].AutoPosition = 0;
View[v0].PositionX = 1000;
View[v0].PositionY = 250;
View[v0].Width = 20;
View[v0].Height = 500;
View[v0].Axes = 2;
View[v0].AxesFormatY = "%.2g";
View[v0].AxesTicsX = 2;
View[v0].AxesTicsY = 5;
// View[v0].ColormapNumber = 5;
// View[v0].ColormapInvert = 1;
// View[v0].ColormapSwap = 1;
Include "map.map";

// Adapt visualization
View[v0].AdaptVisualizationGrid = 1;
View[v0].MaxRecursionLevel = 1;
View[v0].TargetError = -1;

T = t*Dtout;

// Add a title
Plugin(Annotate).Text = Sprintf("t = %5.2f",T);
Plugin(Annotate).X = 1.e5; // By convention, for window coordinates a value greater than 99999 represents the center.
Plugin(Annotate).Y = 70 ;
Plugin(Annotate).FontSize = 24;
Plugin(Annotate).Align = "Center" ;
Plugin(Annotate).View = 0 ;
Plugin(Annotate).Run ;
View[1].Color.Text2D = {255,255,255};

Sleep 0.0001;
Draw;
Print.Background = 1;
Print Sprintf("./pics_rho/kh%03g.jpg", num);

Delete View[1];
Delete View[0];
 
