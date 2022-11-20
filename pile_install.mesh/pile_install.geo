SetFactory("OpenCASCADE");
ancho = 30 ;
prof =-10;
soil1 = -4.7;
soil2= -2.7;
soil3= -11.1;
soil4= -7.7;
soil5= -23.8;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};
//Rectangle(2) = {0, 0, 0, ancho, soil1, 0};

//Rectangle(3) = {0, soil1, 0, ancho, soil2, 0};

//Rectangle(4) = {0, soil1+soil2, 0, ancho, soil3, 0};

//Rectangle(5) = {0, soil1+soil2+soil3, 0, ancho, soil4, 0};

//Rectangle(6) = {0, soil1+soil2+soil3+soil4, 0, ancho, soil5, 0};

//BooleanFragments{ Surface{1}; Delete; }{ Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Delete; }
Physical Line("disp",1) = {4};

Physical Line("far",5) = {2};
Physical Line("under",3) = {3};
Physical Line("level",2) = {1};
Physical Surface("soil1",1)={1};
//Physical Surface("soil2",2)={3};
//Physical Surface("soil3",3)={4};
//Physical Surface("soil4",4)={5};
//Physical Surface("soil5",5)={6};
Transfinite Line { 3} = 70 Using Progression 0.95;
Transfinite Line { 1} = 70 Using Progression 1/0.95;
//+
Transfinite Line{2} =15 Using Progression 1;
Transfinite Line{4} = 120 Using Progression 1;
//+
//Transfinite Surface {1};
//+
//Transfinite Surface {3};
//+
//Transfinite Surface {4};
//+
//Transfinite Surface {5};
//+
//Transfinite Surface {6};
Mesh 2 ;
Mesh.MshFileVersion = 2.2;
RefineMesh;
Save StrCat(StrPrefix(General.FileName), ".msh");