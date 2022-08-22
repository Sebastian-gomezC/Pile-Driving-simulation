// Gmsh project created on Wed Jul 20 11:30:21 2022
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 0.0336, 0.0773, 0};
Physical Line("load",1) = {3};
Physical Line("far",2) = {1};
Physical Surface("soil",1)={1};

