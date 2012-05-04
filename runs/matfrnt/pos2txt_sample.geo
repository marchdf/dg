// Used by pos2txt.sh
// Converts the .pos field data to a ascii .txt file
//
General.Color.Background = White ;
General.Color.Foreground = Black ;
General.BackgroundGradient = 0;
General.SmallAxes = 0;
General.Trackball = 0 ;

Merge Sprintf("../mshORDER/my_line_ORDER_dxDX.msh");
Merge Sprintf("./pORDER/dxDX/FLUX/MODEL/FIELD.pos");

View[0].Axes = 2;
View[0].Color.Axes = Black;
View[0].Type = 2;
View[0].IntervalsType= 2 ; 
View[0].RangeType = 2;
View[0].CustomMax = 1; 
View[0].CustomMin = 0; 
View[0].LineWidth=2;
View[0].AutoPosition = 0;
View[0].Width = 400;
View[0].Height = 260;
View[0].AxesTicsX = 11;
View[0].AxesFormatX = "%.1f";
View[0].AxesTicsY = 10;
View[0].SaturateValues = 1;
View[0].AdaptVisualizationGrid = 1;
View[0].MaxRecursionLevel = ORDER;
View[0].TargetError=0;

Save View[0] Sprintf("./pORDER/dxDX/FLUX/MODEL/FIELD.txt", pdir, dxdir);
Delete View[0];

Exit;
