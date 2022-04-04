
SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0, 0.01};
//+
Point(2) = {1, 0, 0, 0.01};
//+
Line(1) = {1, 2};

Physical Point("left",1) = {1};
Physical Point("right",2) = {2};
Physical Line("Domain",10) = {1};
Mesh 1;
RefineMesh;
RefineMesh;
RefineMesh;
Mesh.MshFileVersion = 2.2;
Save StrCat(StrPrefix(General.FileName), ".msh");
