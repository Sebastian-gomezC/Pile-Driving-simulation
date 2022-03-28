
SetFactory("OpenCASCADE");
ancho = 4;
prof =-1;

Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Physical Curve("disp",1) = {4,2};
Physical Curve("level",2) = {1};
Physical Curve("far",5) = {3};
Physical Surface("soil1",1)={1};

Mesh 2 ;
RefineMesh;
RefineMesh;
Mesh.MshFileVersion = 2.2;
Save StrCat(StrPrefix(General.FileName), ".msh");
