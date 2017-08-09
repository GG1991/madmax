lc=1.5;
lx=1.5;
ly=0.5;
lz=0.5;
r = 0.4;
Nx_cen = 10;
Nx_bor = 10;
Ny_cen = 10;
Ny_bor = 10;
N_cir  = 10;
Nz = 4;
sin45=0.707106781;

Point(1) = {0.0      ,0.0   ,0.0,0.1};
Point(2) = {+lx      ,0.0   ,0.0,0.1};
Point(3) = {+lx      ,+ly   ,0.0,0.1};
Point(4) = {0.0      ,+ly   ,0.0,0.1};
Point(5) = {+ly      ,0.0   ,0.0,0.1};

Circle(1) = {4, 1, 5};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 5};

Transfinite Line {1, 2, 3, 4} = 10 Using Progression 1;

Line Loop(5) = {2, 3, 4, -1};
Plane Surface(6) = {5};
Transfinite Surface {6};
Recombine Surface {6};
