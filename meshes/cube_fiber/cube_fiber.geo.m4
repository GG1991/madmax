//
// Basic 2D figure with hole and a symmetry
//

Geometry.CopyMeshingMethod = 1;

lc=1.5;
lx=1.5;
ly=1.5;
lz=3.0;
r = 0.4;
Nint = Nint_m4;
Next = Nint*2;
Nc = Nint;
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
Line(5) = {3, 4};
Line(6) = {3, 5};

Transfinite Line {1, 3} = Nc Using Progression 1;
Transfinite Line {4, 2} = Next Using Progression 1;
Transfinite Line {5, 6} = Nint Using Progression 1;

loop_m = newl;
Line Loop(loop_m) = {4, 1, 2, 3};

//
// Surface of the matrix
//
surf_m = news;
Plane Surface(surf_m) = {loop_m};
Transfinite Surface {surf_m};
Recombine Surface {surf_m};

//
// Surface of the fiber
//
loop_f = newl;
Line Loop(loop_f) = {6, 1, -5};

surf_f = news;
Plane Surface(surf_f) = {loop_f};
Transfinite Surface {surf_f};
Recombine Surface {surf_f};

//
// First Symmetry
//
sym_1[]=Symmetry {1, -1, 0, 0} {
  Duplicata { Surface{surf_m,surf_f}; }
};

//
// Second Symmetry
//
sym_2[]=Symmetry {0, -1, 0, ly} {
  Duplicata { Surface{surf_m,surf_f,sym_1[0],sym_1[1]}; }
};

//
// Third Symmetry
//
sym_3[]=Symmetry {-1, 0, 0, lx} {
  Duplicata { Surface{surf_m,surf_f,sym_1[0],sym_1[1],sym_2[0],sym_2[1],sym_2[2],sym_2[3]}; }
};

//
// Volume Extrusion
//
vol[]=Extrude {0, 0, lz} {
  Surface{surf_m,surf_f,sym_1[0],sym_1[1],sym_2[0],sym_2[1],sym_2[2],sym_2[3]
    ,sym_3[0],sym_3[1],sym_3[2],sym_3[3],sym_3[4],sym_3[5],sym_3[6],sym_3[7]
  };
  Layers{Nz};
  Recombine;
};

//
// Physical Entities
//
Physical Volume("MATRIX") = {7, 5, 13, 15, 3, 1, 9, 11};
Physical Volume("FIBER") = {8, 6, 4, 2, 10, 12, 16, 14};
Physical Surface("X0") = {127, 205};
Physical Surface("X1") = {361, 283};
Physical Surface("Y0") = {88, 244};
Physical Surface("Y1") = {166, 322};
Physical Surface("Z0") = {20, 53, 62, 44, 35, 8, 11, 29, 34, 25, 58, 67, 49, 40, 10, 16};
Physical Surface("Z1") = {323, 362, 284, 245, 89, 128, 206, 167, 184, 223, 145, 106, 262, 301, 379, 340};
