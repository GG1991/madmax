//
// Basic 2D figure with matrix + fiber
// Non-structured
//

Geometry.CopyMeshingMethod = 1;

lx=3;
ly=3;
dev_x = -0.7;
dev_y = +0.7;
cen_x = 1.5+dev_x;
cen_y = 1.5+dev_y;

r = 0.4;
Nc = 5;
Next = 20;
Nint = 5;

Point(1) = {0.0           ,0.0          ,0.0,0.1};
Point(2) = {+lx           ,0.0          ,0.0,0.1};
Point(3) = {+lx           ,+ly          ,0.0,0.1};
Point(4) = {0.0           ,+ly          ,0.0,0.1};

Point(5) = {cen_x         ,cen_y        ,0.0,0.1};
Point(6) = {cen_x+r       ,cen_y+0      ,0.0,0.1};
Point(7) = {cen_x+0       ,cen_y+r      ,0.0,0.1};
Point(8) = {cen_x-r       ,cen_y+0      ,0.0,0.1};
Point(9) = {cen_x+0       ,cen_y-r      ,0.0,0.1};

Circle(1) = {6, 5, 7};
Circle(2) = {7, 5, 8};
Circle(3) = {8, 5, 9};
Circle(4) = {9, 5, 6};

Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Transfinite Line {8, 7, 6, 5} = 20 Using Progression 1;
Transfinite Line {3, 2, 1, 4} = 5 Using Progression 1;

Line Loop(9) = {3, 4, 1, 2};
Plane Surface(10) = {9};
Line Loop(11) = {8, 5, 6, 7};
Plane Surface(12) = {9, 11};

Physical Surface("FIBER_0")  = {10};
Physical Surface("MATRIX_0") = {12};

trans_1[] = Translate {3.0, 0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_1")  = {trans_1[0]};
Physical Surface("MATRIX_1") = {trans_1[1]};


trans_2[] = Translate {6.0, 0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_2")  = {trans_2[0]};
Physical Surface("MATRIX_2") = {trans_2[1]};


trans_3[] = Translate {9.0, 0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_3")  = {trans_3[0]};
Physical Surface("MATRIX_3") = {trans_3[1]};


trans_4[] = Translate {0, 3.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_4")  = {trans_4[0]};
Physical Surface("MATRIX_4") = {trans_4[1]};


trans_5[] = Translate {3.0, 3.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_5")  = {trans_5[0]};
Physical Surface("MATRIX_5") = {trans_5[1]};


trans_6[] = Translate {6.0, 3.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_6")  = {trans_6[0]};
Physical Surface("MATRIX_6") = {trans_6[1]};


trans_7[] = Translate {9.0, 3.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_7")  = {trans_7[0]};
Physical Surface("MATRIX_7") = {trans_7[1]};


trans_8[] = Translate {0, 6.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_8")  = {trans_8[0]};
Physical Surface("MATRIX_8") = {trans_8[1]};


trans_9[] = Translate {3.0, 6.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_9")  = {trans_9[0]};
Physical Surface("MATRIX_9") = {trans_9[1]};


trans_10[] = Translate {6.0, 6.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_10")  = {trans_10[0]};
Physical Surface("MATRIX_10") = {trans_10[1]};


trans_11[] = Translate {9.0, 6.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_11")  = {trans_11[0]};
Physical Surface("MATRIX_11") = {trans_11[1]};


trans_12[] = Translate {0, 9.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_12")  = {trans_12[0]};
Physical Surface("MATRIX_12") = {trans_12[1]};


trans_13[] = Translate {3.0, 9.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_13")  = {trans_13[0]};
Physical Surface("MATRIX_13") = {trans_13[1]};


trans_14[] = Translate {6.0, 9.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_14")  = {trans_14[0]};
Physical Surface("MATRIX_14") = {trans_14[1]};


trans_15[] = Translate {9.0, 9.0, 0} {
  Duplicata { Surface{10, 12}; }
};

Physical Surface("FIBER_15")  = {trans_15[0]};
Physical Surface("MATRIX_15") = {trans_15[1]};
        
Physical Line("Y0") = {5, 25, 38, 51};
Physical Line("Y1") = {162, 176, 188, 200};
Physical Line("X0") = {8, 65, 115, 165};
Physical Line("X1") = {50, 101, 151, 201};
