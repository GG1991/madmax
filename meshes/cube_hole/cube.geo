lc=1.5;
lx=0.5;
ly=0.5;
r = 0.25;
Nx_cen = 10;
Nx_bor = 10;
Ny_cen = 10;
Ny_bor = 10;
N_cir  = 10;

Point(1) = {-lx   ,-ly   ,-lc,0.1};
Point(2) = {+lx   ,-ly   ,-lc,0.1};
Point(3) = {+lx   ,+ly   ,-lc,0.1};
Point(4) = {-lx   , ly   ,-lc,0.1};

Point(5) = {-lx+r ,-ly   ,-lc,0.1};
Point(6) = {-lx+r ,-ly+r ,-lc,0.1};
Point(7) = {-lx   ,-ly+r ,-lc,0.1};

Point(8) = {+lx-r ,+ly   ,-lc,0.1};
Point(9) = {+lx-r ,+ly-r ,-lc,0.1};
Point(10)= {+lx   ,+ly-r ,-lc,0.1};

Point(11) = {-lx+r ,+ly   ,-lc,0.1};
Point(12) = {-lx+r ,+ly-r ,-lc,0.1};
Point(13) = {-lx   ,+ly-r ,-lc,0.1};
Point(14) = {+lx-r ,-ly+r ,-lc,0.1};
Point(15) = {+lx-r ,-ly   ,-lc,0.1};
Point(16) = {+lx   ,-ly+r ,-lc,0.1};

// arcs
Circle(1) = {7, 1, 5};
Circle(2) = {8, 3, 10};
Transfinite Line {1, 2} = N_cir Using Progression 1;

//lines horizontal border
Line(3) = {7, 6};
Line(4) = {9, 10};
Line(5) = {13, 12};
Line(6) = {4, 11};
Line(7) = {14, 16};
Line(8) = {15, 2};
Transfinite Line {3, 5, 6, 4, 7, 8} = Nx_bor Using Progression 1;

//lines vertical border
Line(9) = {4, 13};
Line(10) = {11, 12};
Line(11) = {8, 9};
Line(12) = {6, 5};
Line(13) = {14, 15};
Line(14) = {16, 2};
Transfinite Line {9, 10, 11, 12, 13, 14} = Ny_bor Using Progression 1;

//lines horizontal center
Line(15) = {12, 9};
Line(16) = {11, 8};
Line(17) = {6, 14};
Line(18) = {5, 15};
Transfinite Line {15,16,17,18} = Nx_cen Using Progression 1;

//lines vertical center
Line(19) = {13, 7};
Line(20) = {12, 6};
Line(21) = {9, 14};
Line(22) = {10, 16};
Transfinite Line {19,20,21,22} = Ny_cen Using Progression 1;

// surfaces
Line Loop(23) = {19, 3, -20, -5};
Plane Surface(24) = {23};
Transfinite Surface {24};

Line Loop(25) = {21, 7, -22, -4};
Plane Surface(26) = {25};
Transfinite Surface {26};

Line Loop(27) = {17, 13, -18, -12};
Plane Surface(28) = {27};
Transfinite Surface {28};

Line Loop(29) = {10, 15, -11, -16};
Plane Surface(30) = {29};
Transfinite Surface {30};

Line Loop(31) = {9, 5, -10, -6};
Plane Surface(32) = {31};
Transfinite Surface {32};

Line Loop(33) = {13, 8, -14, -7};
Plane Surface(34) = {33};
Transfinite Surface {34};

Line Loop(35) = {17, -21, -15, 20};
Plane Surface(36) = {35};
Transfinite Surface {36};

Line Loop(37) = {3, 12, -1};
Plane Surface(38) = {37};
Transfinite Surface {38};

Line Loop(39) = {11, 4, -2};
Plane Surface(40) = {39};
Transfinite Surface {40};
