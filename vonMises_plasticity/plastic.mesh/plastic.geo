SetFactory("OpenCASCADE");
ancho = 10;
prof =-3;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Rectangle(2) = {ancho/2-0.2, 0, 0,0.4, -0.4, 0};
Rectangle(3) = {ancho/2-0.5, -0.4, 0,1, -0.2, 0};
BooleanDifference{ Surface{1}; Delete; }{ Surface{3}; Surface{2}; Delete; }


Physical Line("disp",1) = {13,15,12,10,8,6};
Physical Line("load",2) = {11};
Physical Line("level",3) = {19,16};
Physical Line("far",5) = {14};
Physical Surface("soil1",1)={1};
Transfinite Curve{8,18,12,10,17,6,16,19} = 25 ;
Transfinite Curve{19} = 60 Using Progression 0.98;
Transfinite Curve{16} = 60 Using Progression 1/0.98;
Transfinite Curve{11} = 50 ;
Transfinite Curve {14} = 80 Using Bump 1.5;
Transfinite Curve{15,13} = 20 ;
Mesh 2 ;
Mesh.MshFileVersion = 2.2;
