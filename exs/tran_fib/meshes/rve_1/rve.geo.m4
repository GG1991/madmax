//
// Basic 2D figure with hole and a symmetry
//

Geometry.CopyMeshingMethod = 1;

lc=1.5;
lx=1.5;
ly=1.5;
r = 0.4;
Nint = n_m4;
Next = Nint*2;
Nc   = Nint;
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

Physical Surface("FIBER")  = {34, 25, 10, 40, 16, 58, 49, 67};
Physical Surface("MATRIX") = {29, 20, 53, 62, 44, 35, 8, 11};
Physical Line("X0") = {15, 33};
Physical Line("X1") = {48, 66};
Physical Line("Y0") = {3, 39};
Physical Line("Y1") = {24, 57};
