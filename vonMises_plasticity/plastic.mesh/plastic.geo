SetFactory("OpenCASCADE");
ancho = 5;
prof =3;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Rectangle(2) = {0, prof, 0,ancho/2-0.3, -0.1, 0};
Rectangle(3) = {ancho/2+0.3, prof, 0,ancho/2-0.3, -0.1, 0};
BooleanDifference{ Surface{1}; Delete; }{ Surface{3}; Surface{2}; Delete; }
Rectangle(4) = {ancho/2-0.3,prof,0,0.6,-0.1,0};
BooleanDifference{ Surface{1}; Delete; }{ Surface{4}; Delete; }
Physical Line("disp",1) = {6,4};
Physical Line("load",2) = {2};
//Physical Line("level",3) = {1,5,6};
Physical Line("far",5) = {5};
Physical Surface("soil",1)={1};
Transfinite Curve {2} = 80 ;
Transfinite Curve {4} = 40 Using Progression 0.98 ;
Transfinite Curve {6} = 40 Using Progression 1.02 ;
Transfinite Curve {5} = 20;
Transfinite Curve {1,3} = 60;
Mesh 2 ;
Mesh.MshFileVersion = 2.2;
