
SetFactory("OpenCASCADE");
ancho = 30 ;
prof =-50;
soil1 = -4.7;
soil2= -2.7;
soil3= -11.1;
soil4= -7.7;
soil5= -23.8;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};


Rectangle(2) = {0, 0, 0, ancho, soil1, 0};

Rectangle(3) = {0, soil1, 0, ancho, soil2, 0};

Rectangle(4) = {0, soil1+soil2, 0, ancho, soil3, 0};

Rectangle(5) = {0, soil1+soil2+soil3, 0, ancho, soil4, 0};

Rectangle(6) = {0, soil1+soil2+soil3+soil4, 0, ancho, soil5, 0};

BooleanFragments{ Surface{1}; Delete; }{ Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Delete; }
Physical Curve("disp",1) = {1,5,8,11,14};
Physical Curve("level",2) = {4};
Physical Curve("far",5) = {3,7,10,13,15,16};
Physical Surface("soil1",1)={2};
Physical Surface("soil2",2)={3};
Physical Surface("soil3",3)={4};
Physical Surface("soil4",4)={5};
Physical Surface("soil5",5)={6};
Characteristic Length {4,3,6,8,10,12} = 3;
Characteristic Length {2,1,5,7,9,11} = 0.1;
Mesh 2 ;

Mesh.MshFileVersion = 2.2;
Save StrCat(StrPrefix(General.FileName), ".msh");
