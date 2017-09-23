//
// Basic 2D figure with hole and a symmetry
//

Geometry.CopyMeshingMethod = 1;

lx =30;
ly =lx;
Nx = n_m4;
Ny = Nx;

Point(1) = {0.0  ,0.0 ,0.0,0.1};
Point(2) = {+lx  ,0.0 ,0.0,0.1};
Point(3) = {+lx  ,+ly ,0.0,0.1};
Point(4) = {0.0  ,+ly ,0.0,0.1};

Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};

Transfinite Line {3, 1} = Nx Using Progression 1;
Transfinite Line {4, 2} = Ny Using Progression 1;
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Transfinite Surface {6};
Recombine Surface {6};

Physical Surface("MICRO") = {6};
Physical Line("Y0") = {3};
Physical Line("Y1") = {1};
Physical Line("X0") = {4};
Physical Line("X1") = {2};
