//
// Basic 2D figure with hole and a symmetry
//

Geometry.CopyMeshingMethod = 1;

lc=1.5;
lx=1.5;
ly=1.5;
lz=0.5;
r = 0.4;
Nc = 2;
Nx = 2;
Nz = 1;
sin45=0.707106781;

Point(1) = {0.0           ,0.0          ,0.0,0.1};
Point(2) = {+lx           ,0.0          ,0.0,0.1};
Point(3) = {+lx           ,+ly          ,0.0,0.1};
Point(4) = {+lx           ,+ly-r        ,0.0,0.1};
Point(5) = {+lx-r*sin45   ,+ly-r*sin45  ,0.0,0.1};

Circle(1) = {5, 3, 4};
Line(2) = {4, 2};
Line(3) = {2, 1};
Line(4) = {1, 5};
Transfinite Line {1, 3} = Nc Using Progression 1;
Transfinite Line {4, 2} = Nx Using Progression 1;

lloop_1 = newl;
Line Loop(lloop_1) = {4, 1, 2, 3};

surf_1 = news;
Plane Surface(surf_1) = {5};
Transfinite Surface {surf_1};
Recombine Surface {surf_1};

surf_2[]=Symmetry {1, -1, 0, 0} {
  Duplicata { Surface{surf_1}; }
};

surf_3[]=Symmetry {0, -1, 0, ly} {
  Duplicata { Surface{surf_1,surf_2[0]}; }
};

surf_4[]=Symmetry {-1, 0, 0, lx} {
  Duplicata { Surface{surf_1,surf_2[0],surf_3[0],surf_3[1]}; }
};

//
// Volume Extrusion
//
vol[]=Extrude {0, 0, 1} {
  Surface{surf_1,surf_2[0],surf_3[0],surf_3[1],surf_4[0],surf_4[1],surf_4[2],surf_4[3]};
  Layers{Nz};
  Recombine;
};

//
// Physical Entities
//
Physical Volume("IRON") = {3, 4, 2, 1, 5, 6, 8, 7};
Physical Surface("X0") = {128, 84};
Physical Surface("X1") = {216, 172};
Physical Surface("Y0") = {62, 150};
Physical Surface("Y1") = {106, 194};
Physical Surface("Z0") = {12, 17, 7, 6, 22, 27, 37, 32};
Physical Surface("Z1") = {107, 129, 85, 63, 151, 173, 217, 195};
