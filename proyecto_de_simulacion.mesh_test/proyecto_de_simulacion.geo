
    SetFactory("OpenCASCADE");
    h=5;
    l=15;
    Point(1) = {0,0,0};
    Point(2) = {0,20,0};
    Point(3) = {15-l/2,20,0};
    Point(4) = {15-l/2,20-h,0};
    Point(5) = {15+l/2,20-h,0};
    Point(6) = {15+l/2-1,20,0};
    Point(7) = {30,20,0};
    Point(8) = {30,0,0};
    
    Line(1) = {1, 2};
    
    Line(2) = {2, 3};
    
    Line(3) = {3, 4};
    //+
    Line(4) = {4, 5};
    //+
    Line(5) = {5, 6};
    //+
    Line(6) = {6, 7};
    //+
    Line(7) = {7, 8};
    //+
    Line(8) = {8, 1};
    //+
    Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
    //+
    Plane Surface(1) = {1};
    Physical Curve("a",1)={2};
    Physical Curve("bordes",3)={1,7};
    Physical Curve("neuman0",4)={3,4,5,8};
    Physical Curve("c",6)={6};
    Physical Surface("superficie",5)={1};
    
    Mesh 2;
    Coherence Mesh;
    
    RefineMesh;
    Mesh.MshFileVersion = 2.2;
    Save StrCat(StrPrefix(General.FileName), ".msh");
    
