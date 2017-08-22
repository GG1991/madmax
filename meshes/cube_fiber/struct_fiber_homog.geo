//
// Basic 2D figure with hole and a symmetry
//

Geometry.CopyMeshingMethod = 1;

lc=1.5;
lx=1.5;
ly=1.5;
lz=lx*2*3;
Nx = 4;
Ny = 4;
Nz = 1;

Point(1) = {0.0      ,0.0     ,0.0,0.1};
Point(2) = {+lx*2*3  ,0.0     ,0.0,0.1};
Point(3) = {+lx*2*3  ,+ly*2*3 ,0.0,0.1};
Point(4) = {0.0      ,+ly*2*3 ,0.0,0.1};

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

Extrude {0, 0, lz} {
  Surface{6};
  Layers{Nz};
  Recombine;
}

Physical Volume("MICRO") = {1};
Physical Surface("X0") = {15};
Physical Surface("X1") = {23};
Physical Surface("Y0") = {27};
Physical Surface("Y1") = {19};
Physical Surface("Z0") = {6};
Physical Surface("Z1") = {28};
