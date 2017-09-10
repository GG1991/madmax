lc=1.5;
N = 100;

Point(1) = {-lc,-lc,-lc,0.1};
Point(2) = {lc,-lc,-lc,0.1};
Point(3) = {lc,lc,-lc,0.1};
Point(4) = {-lc, lc,-lc,0.1};

Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,1};
Line(4) = {1,2};

Transfinite Line {1, 3} = N Using Progression 1;
Transfinite Line {2, 4} = N Using Progression 1;
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
Transfinite Surface {6};
Recombine Surface {6};

Physical Point("P000") = {1};
Physical Point("P100") = {2};
Physical Point("P010") = {4};

Physical Surface("MATRIX") = {6};
Physical Line("X0") = {3};
Physical Line("X1") = {1};
Physical Line("Y0") = {4};
Physical Line("Y1") = {2};
