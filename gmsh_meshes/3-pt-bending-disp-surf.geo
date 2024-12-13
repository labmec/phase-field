// Gmsh project created on Thu Feb  1 11:32:44 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {3.75, 0, 0, 1.0};
Point(3) = {4.25, 0, 0, 1.0};
Point(4) = {8, 0, 0, 1.0};
Point(5) = {8, 2, 0, 1.0};
Point(6) = {4.25, 2, 0, 1.0};
Point(7) = {3.75, 2, 0, 1.0};
Point(8) = {0, 2, 0, 1.0};
Point(9) = {3.95, 0, 0, 1.0};
Point(10) = {4.05, 0, 0, 1.0};
Point(11) = {4.0, 0.2, 0, 1.0};
Point(12) = {4.0, 2.0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 9};
//+
Line(3) = {9, 11};
//+
Line(4) = {11, 10};
//+
Line(5) = {10, 3};
//+
Line(6) = {3, 4};
//+
Line(7) = {4, 5};
//+
Line(8) = {5, 6};
//+
Line(9) = {6, 12};
//+
Line(10) = {12, 7};
//+
Line(11) = {7, 8};
//+
Line(12) = {8, 1};
//+
Line(13) = {7, 2};
//+
Line(14) = {6, 3};
//+
Line(15) = {12, 11};
//+
Curve Loop(1) = {1, -13, 11, 12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 7, 8, 14};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {4, 5, -14, 9, 15};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, -15, 10, 13, 2};
//+
Plane Surface(4) = {4};
//+
Physical Point("ptfixedxy", 16) = {1};
//+
Physical Point("ptfixedy", 17) = {4};
//+
Physical Surface("dom", 19) = {1, 2};
//+
Physical Curve("dispy", 21) = {10, 9};
//+
Physical Surface("domfrac", 20) = {4, 3};
//+
Recombine Surface {1, 2};
Recombine Surface {3, 4};

//+ Top and bottom of left side of domain
Transfinite Curve {1, 11} = 15 Using Progression 1;

//+ Left and right of left side of domain
Transfinite Curve {12} = 8 Using Progression 1;
Transfinite Curve {13} = 40 Using Progression 1;
//+ Left and right of right side of domain
Transfinite Curve {7} = 8 Using Progression 1;
Transfinite Curve {14} = 40 Using Progression 1;

//+ Top and bottom of right side of domain
Transfinite Curve {6, 8} = 15 Using Progression 1;

//+ Middle where fracture is
Transfinite Curve {15} = 1000 Using Progression 1;
//+ Top in the middle part of domain
Transfinite Curve {10, 9} = 10 Using Progression 1;
//+ Bottom in the middle part of domain
Transfinite Curve {2, 5} = 8 Using Progression 1;
//+ Notch edges
Transfinite Curve {3, 4} = 80 Using Progression 1;

