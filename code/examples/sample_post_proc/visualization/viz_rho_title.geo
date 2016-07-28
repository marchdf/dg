//
// Another example of visualizing the density field. This one was
// specifically designed to generate a pretty picture for a title
// slide.
// 

Merge Sprintf("./rho%010.0f.pos",num);

// Assign values to the views
v0 = PostProcessing.NbViews-1;

// Custom density
View[v0].Name = "Density";
View[v0].ScaleType = 1; 
View[v0].RangeType = 2;
View[v0].CustomMax = 1;
View[v0].CustomMin = 0.33;
View[v0].SaturateValues = 1;
View.ShowScale = 0;
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
// View[v0].Axes = 2;
// View[v0].AxesFormatY = "%.2g";
// View[v0].AxesTicsX = 2;
// View[v0].AxesTicsY = 5;
// View[v0].ColormapNumber = 5;
// View[v0].ColormapInvert = 1;
// View[v0].ColormapSwap = 1;
Include "map.map";

// Adapt visualization
View[v0].AdaptVisualizationGrid = 1;
View[v0].MaxRecursionLevel = 1;
View[v0].TargetError = -1;

// Clip
General.Clip0A = 0;
General.Clip0B = 1;
General.Clip0C = 0;
General.Clip0D = 0.5; 
General.Clip1A = 0;
General.Clip1B =-1;
General.Clip1C = 0;
General.Clip1D = 0.5;

View[v0].Clip=3; // 2^0 + 2^1


// // Alias the views
// AliasWithOptions View[v0];
// AliasWithOptions View[v0];
// AliasWithOptions View[v0];
// AliasWithOptions View[v0];
// AliasWithOptions View[v0];

// // Shift them
// View[v0+1].OffsetX = -1;
// View[v0+2].OffsetX = 1;
// View[v0+3].OffsetX = -2;
// View[v0+4].OffsetX = 2;
// View[v0+5].OffsetX = 3;

// turn off the lighting
View[v0].Light = 0;
// View[v0+1].Light = 0;
// View[v0+2].Light = 0;
// View[v0+3].Light = 0;
// View[v0+4].Light = 0;
// View[v0+5].Light = 0;



Sleep 0.0001;
Draw;
Print.Background = 1;
Print Sprintf("./pics_rho/kh%03g.jpg", num);

Delete View[0];
