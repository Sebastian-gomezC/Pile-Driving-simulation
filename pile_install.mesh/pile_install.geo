

SetFactory("OpenCASCADE");
ancho = 30 ;
prof =-10;
soil1 = -4.7;
soil2= -2.7;
soil3= -11.1;
soil4= -7.7;
soil5= -23.8;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Physical Curve("disp",1) = {4};
Physical Curve("level",2) = {1};
Physical Curve("far",5) = {3,2};
Physical Surface("soil1",1)={1};
Characteristic Length {2,3} = 1;
Characteristic Length {1,4} = 0.08;
Mesh 2 ;

Mesh.MshFileVersion = 2.2;
Save StrCat(StrPrefix(General.FileName), ".msh");
