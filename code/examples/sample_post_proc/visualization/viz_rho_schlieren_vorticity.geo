//
// Visualize density, schlieren and vorticity in a subset of the
// domain.  Loops over all time and saves images to the pics/
// directory (should exist prior to running this script)
//

Merge "./mesh.msh";
Merge "./rho.pos";
Merge "./ux.pos";
Merge "./uy.pos";
Dtout = 1e-4;
tConvert = 1e3;

General.Color.Background = Black ;
General.Color.Foreground = White ;
General.Color.Text = White ;
General.BackgroundGradient = 0;
General.SmallAxes = 0;
General.Trackball = 0 ;
General.RotationX = 0;
General.RotationY = 0;
General.RotationZ = 00;
General.ScaleX = 2.4;
General.ScaleY = 2.4;
General.TranslationX = 0;
General.TranslationY = 0.01;
General.GraphicsWidth = 2400 ;
General.GraphicsHeight = 1280 ;
General.GraphicsFontSize = 18;
General.GraphicsFontSizeTitle = 24;

// Cut the views to the domain size we want
General.Clip0A = 0;
General.Clip0B = 1;
General.Clip0C = 0;
General.Clip0D = 0.5;

Geometry.Points = 0;
Geometry.Clip = 1;

Mesh.Lines = 0;
Mesh.LineWidth = 2;
Mesh.Triangles = 1;
Mesh.SurfaceEdges=0;
Mesh.SurfaceFaces=0;
Mesh.ColorCarousel = 0; //Mesh coloring (0=by element type, 1=by elementary entity, 2=by physical entity, 3=by partition)
Mesh.Clip=1;

PostProcessing.HorizontalScales=0;

// Assign values to the views
v0 = PostProcessing.NbViews-3;
v1 = PostProcessing.NbViews-2;
v2 = PostProcessing.NbViews-1;

// Custom density
View[v0].RangeType = 2;
View[v0].CustomMax = 10;
View[v0].CustomMin = 1.2;
View[v0].SaturateValues = 1;
View[v0].Format = "%5.5g";
View[v0].NbIso = 4;
View[v0].ShowTime = 0;
View[v0].Color.Axes = White;
View[v0].OffsetX = -0.07;
View[v0].AutoPosition = 0;
View[v0].PositionX = 460;
View[v0].PositionY = 150;
View[v0].Width = 20;
View[v0].Height = 430;
View[v0].Axes = 1;
View[v0].AxesFormatY = "%.2g";
View[v0].AxesTicsX = 0;
View[v0].AxesTicsY = 3;
Include "map.map";

// Make velocities invisible
View[v1].Visible = 0;
View[v2].Visible = 0;

// View 3: make a vector view of velocities
Plugin(Scal2Vec).NameNewView="U";
Plugin(Scal2Vec).ViewX=v1;
Plugin(Scal2Vec).ViewY=v2;
Plugin(Scal2Vec).Run;

Delete View[v2];
Delete View[v1]; // View 3 is now view 1

// View 2: take the curl of velocity
Plugin(Curl).View = 1;
Plugin(Curl).Run;

Delete View[1]; // View 2 is now view 1

// View 2: take the 3d component of vorticity
Plugin(MathEval).Expression0 = "v2";
Plugin(MathEval).TimeStep = -1;
Plugin(MathEval).View = 1;
Plugin(MathEval).Run;
Delete View[1];
// View 2 is now view 1
View[1].Name = "Vort.";
View[1].RangeType = 2;
View[1].CustomMax = 3000;
View[1].CustomMin =-3000;
View[1].NbIso = 5;
View[1].SaturateValues = 1;
View[1].Format = "%5.5g";
View[1].ShowTime = 0;
View[1].Color.Axes = White;
// View[1].AdaptVisualizationGrid = 1;
// View[1].MaxRecursionLevel = 1;
// View[1].TargetError = -1;
View[1].OffsetX = 0.07;
View[1].AutoPosition = 0;
View[1].PositionX = 990;
View[1].PositionY = 150;
View[1].Width = 20;
View[1].Height = 430;
Include "map1.map";

// View 2: take the gradient of density
Plugin(Gradient).View = v0;
Plugin(Gradient).Run;                  // run the plugin!

// View 3: take the norm of density
Plugin(MathEval).Expression0 = "Sqrt(v0^2+v1^2+v2^2)";
Plugin(MathEval).TimeStep = -1;
Plugin(MathEval).View = 2;
Plugin(MathEval).Run;

Delete View[2];
// View 3 is now view 2

// Get Schlieren view

// Find the Maximum of the norm of the gradient
Max = View[2].Max;

// View 3: Take the norm of the gradient of density
Plugin(MathEval).Expression0 = Sprintf("Exp(-20*v0/%e)",Max);
Plugin(MathEval).TimeStep = -1;
Plugin(MathEval).View = 2;
Plugin(MathEval).Run;

Delete View[2]; // View 3 is now view 2

// Custom the Schlieren view
View[v2].Name = "Sch.";
View[v2].RangeType = 2;
View[v2].CustomMax = 1;
View[v2].CustomMin = .01;
View[v2].OffsetX = 0;
View[v2].SaturateValues = 1;
View[v2].ColormapNumber = 9;
View[v2].ShowScale = 0;

// Cut the views to the domain size we want
View[v0].Clip = 1;
View[v1].Clip = 1;
View[v2].Clip = 1;


// Loop on time to get the Schlieren view
t = 0;
Nt = View[v0].NbTimeStep;

// View[v0].AdaptVisualizationGrid = 1;
// View[v0].MaxRecursionLevel = 1;
// View[v0].TargetError = -1;

For num In {0:Nt-1}

// Plot the density and vorticity
View[v0].TimeStep = t ;
View[v1].TimeStep = t ;
View[v2].TimeStep = t ;

T = t*Dtout*tConvert;

// Add a title
Plugin(Annotate).Text = Sprintf("t = %4.1f [ms]",T);
Plugin(Annotate).X = 1.e5; // By convention, for window coordinates a value greater than 99999 represents the center.
Plugin(Annotate).Y = 70 ;
Plugin(Annotate).FontSize = 24;
Plugin(Annotate).Align = "Center" ;
Plugin(Annotate).View = 0 ;
Plugin(Annotate).Run ;
View[3].Color.Text2D = {255,255,255};

Sleep 0.0001;
Draw;
Print.Background = 1;
Print Sprintf("./pics2/rm%03g.jpg", num);

t = t+1;

Delete View[3];

EndFor

// Exit;
// //System "mencoder mf://pics/*.jpg -mf fps=10:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi";
// //System "ffmpeg -i output.avi -vcodec mpeg4 -b:v 1200k -flags +aic+mv4 output.mp4"; // from http://andrebluehs.net/blog/2012/05/converting-avi-to-mp4-with-ffmpeg/
