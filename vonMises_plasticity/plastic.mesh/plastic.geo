SetFactory("OpenCASCADE");
ancho = 10;
prof =-3;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Rectangle(2) = {ancho/2-0.2, 0, 0,0.4, -0.4, 0};
Rectangle(3) = {ancho/2-0.5, -0.4, 0,1, -0.2, 0};
BooleanDifference{ Surface{1}; Delete; }{ Surface{3}; Surface{2}; Delete; }
Physical Line("disp",1) = {9,11,10,5,7,6};
Physical Line("load",2) = {8};
Physical Line("level",3) = {12,4};
Physical Line("far",5) = {2};
Physical Surface("soil1",1)={1};
Transfinite Curve{9,11,10,5,7,6} = 25 ;
Transfinite Curve{12} = 60 Using Progression 0.98;
Transfinite Curve{4} = 60 Using Progression 1/0.98;
Transfinite Curve{8} = 50 ;
Transfinite Curve {2} = 80 Using Bump 1.5;
Transfinite Curve{1,3} = 10 ;
Mesh 2 ;
Mesh.MshFileVersion = 2.2;
