//
// Basic 2D figure with hole and a symmetry
//

Geometry.CopyMeshingMethod = 1;

lx=3.0;
ly=3.0;
Ny = n_m4 ;
Nx = Ny * 3;

Point(1) = {0.0           ,0.0          ,0.0,0.1};
Point(2) = {+lx           ,0.0          ,0.0,0.1};
Point(3) = {+lx           ,+ly          ,0.0,0.1};
Point(4) = {0.0           ,+ly          ,0.0,0.1};

Point(5) = {0.0           ,+ly*1.0/3.0  ,0.0,0.1};
Point(6) = {0.0           ,+ly*2.0/3.0  ,0.0,0.1};
Point(7) = {+lx           ,+ly*1.0/3.0  ,0.0,0.1};
Point(8) = {+lx           ,+ly*2.0/3.0  ,0.0,0.1};


Line(1) = {4, 6};
Line(2) = {6, 5};
Line(3) = {5, 1};
Line(4) = {1, 2};
Line(5) = {2, 7};
Line(6) = {7, 8};
Line(7) = {8, 3};
Line(8) = {3, 4};
Line(9) = {6, 8};
Line(10) = {5, 7};

Line Loop(11) = {1, 9, 7, 8};
Plane Surface(12) = {11};
Line Loop(13) = {2, 10, 6, -9};
Plane Surface(14) = {13};
Line Loop(15) = {10, -5, -4, -3};
Plane Surface(16) = {15};
Transfinite Line {2, 3, 1, 7, 6, 5} = Ny Using Progression 1;
Transfinite Line {4, 10, 9, 8}      = Nx Using Progression 1;
Transfinite Surface {16};
Transfinite Surface {14};
Transfinite Surface {12};

Recombine Surface {16, 14, 12};

Physical Line("Y0") = {4};
Physical Line("Y1") = {8};
Physical Line("X0") = {3, 2, 1};
Physical Line("X1") = {5, 6, 7};
Physical Surface("FIBER")  = {14};
Physical Surface("MATRIX") = {12, 16};
