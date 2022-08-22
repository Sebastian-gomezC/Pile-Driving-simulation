SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 0.0336, 0.0773, 0};
//+
Physical Curve("load", 1) = {3};
//+
Physical Curve("far", 2) = {1};
//+
Physical Curve("borders", 3) = {4, 2};
//+
}Physical Surface("soil", 1) = {1};
Transfinite Curve {3, 1} = 25 Using Progression 1;
//+
Transfinite Curve {4, 2} = 50 Using Progression 1;
//+
Transfinite Surface {1};
Mesh 2 ;
Mesh.MshFileVersion = 2.2;
