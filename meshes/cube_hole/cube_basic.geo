lc=1.5;
lx=1.5;
ly=1.5;
lz=0.5;
r = 0.4;
Nx_cen = 10;
Nx_bor = 10;
Ny_cen = 10;
Ny_bor = 10;
N_cir  = 10;
Nz = 4;
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
Transfinite Line {1, 3} = 10 Using Progression 1;
Transfinite Line {4, 2} = 10 Using Progression 1;

Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Transfinite Surface {6};
Recombine Surface {6};
