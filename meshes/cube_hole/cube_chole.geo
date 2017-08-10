lc=0.5;
lx=0.5;
ly=0.5;
lh=0.2;
N=1;
NH=3*N;
NR=9*N;

Point(1) = {-lx   ,-ly   ,-lc,0.1};
Point(2) = {+lx   ,-ly   ,-lc,0.1};
Point(3) = {+lx   ,+ly   ,-lc,0.1};
Point(4) = {-lx   , ly   ,-lc,0.1};

Point(5) = {+lx-lh   ,+ly    ,-lc,0.1};
Point(6) = {+lx-lh   ,+ly-lh ,-lc,0.1};
Point(7) = {+lx-lh   ,-ly    ,-lc,0.1};

Point(8) = {+lx      ,+ly-lh ,-lc,0.1};
Point(9) = {-lx      ,+ly-lh ,-lc,0.1};

Line(1) = {9, 6};
Line(2) = {6, 7};
Line(3) = {7, 1};
Line(4) = {1, 9};
Line(5) = {4, 5};
Line(6) = {3, 8};
Line(7) = {5, 6};
Line(8) = {4, 9};
Line(9) = {5, 3};
Line(10) = {8, 6};
Line(11) = {8, 2};
Line(12) = {2, 7};
Transfinite Line {6, 7, 8, 9, 10, 12} = NH Using Progression 1;
Transfinite Line {5, 1, 2, 11, 3, 4} = NR Using Progression 1;
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};
Line Loop(15) = {1, -7, -5, 8};
Plane Surface(16) = {15};
Line Loop(17) = {10, 2, -12, -11};
Plane Surface(18) = {17};
Transfinite Surface {14};
Transfinite Surface {16};
Transfinite Surface {18};
Recombine Surface {14};
Recombine Surface {16};
Recombine Surface {18};

vol[]=Extrude {0, 0, 2*lc} {
  Surface{14, 16, 18};
  Layers{10*N};
  Recombine;
};

Physical Volume("IRON") = {1, 2, 3};
Physical Surface("X0") = {39, 61};
Physical Surface("X1") = {83};
Physical Surface("Y0") = {35, 79};
Physical Surface("Y1") = {57};
Physical Surface("Z0") = {16, 14, 18};
Physical Surface("Z1") = {62, 40, 84};
