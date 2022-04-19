
SetFactory("OpenCASCADE");
ancho = 0.5;
prof =-0.22;

Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Physical Line("disp",1) = {4,2};
Physical Line("level",2) = {1};
Physical Line("far",5) = {3};
Physical Surface("soil1",1)={1};

Mesh 2 ;
RefineMesh;
RefineMesh;
RefineMesh;
Mesh.MshFileVersion = 2.2;
Save StrCat(StrPrefix(General.FileName), ".msh");
