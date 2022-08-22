
SetFactory("OpenCASCADE");
ancho = 0.05;
prof =-0.1;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Physical Line("disp",1) = {4,2};
Physical Line("level",2) = {1};
Physical Line("far",5) = {3};
Physical Surface("soil1",1)={1};
Transfinite Curve{1} = 80 ;
Transfinite Curve {2} = 40 Using Progression 1.05;
Transfinite Curve {4} = 40 Using Progression 1/1.05;
Mesh 2 ;
Mesh.MshFileVersion = 2.2;
Save StrCat(StrPrefix(General.FileName), ".msh");
