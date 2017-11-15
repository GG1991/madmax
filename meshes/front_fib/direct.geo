//
// Basic 2D figure with hole and a symmetry
//

Geometry.CopyMeshingMethod = 1;

lc=1.5;
lx=1.5;
ly=1.5;
lz=lx*2*3;
r = 0.4;
Nc = 3;
Next = 10;
Nint = 3;
Nz = 1;
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

Physical Surface("FIBER_0") ={surf_f,sym_1[1],sym_2[1],sym_2[3],sym_3[1],sym_3[3],sym_3[5],sym_3[7]};
Physical Surface("MATRIX_0")={surf_m,sym_1[0],sym_2[0],sym_2[2],sym_3[0],sym_3[2],sym_3[4],sym_3[6]};


tras_1[]=Translate {lx*2*0, ly*2*1, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_1") ={tras_1[4],tras_1[5],tras_1[6],tras_1[7],tras_1[8],tras_1[9],tras_1[10],tras_1[11]};
Physical Surface("MATRIX_1") ={tras_1[0],tras_1[1],tras_1[2],tras_1[3],tras_1[12],tras_1[13],tras_1[14],tras_1[15]};


tras_2[]=Translate {lx*2*0, ly*2*2, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_2") ={tras_2[4],tras_2[5],tras_2[6],tras_2[7],tras_2[8],tras_2[9],tras_2[10],tras_2[11]};
Physical Surface("MATRIX_2") ={tras_2[0],tras_2[1],tras_2[2],tras_2[3],tras_2[12],tras_2[13],tras_2[14],tras_2[15]};


tras_3[]=Translate {lx*2*0, ly*2*3, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_3") ={tras_3[4],tras_3[5],tras_3[6],tras_3[7],tras_3[8],tras_3[9],tras_3[10],tras_3[11]};
Physical Surface("MATRIX_3") ={tras_3[0],tras_3[1],tras_3[2],tras_3[3],tras_3[12],tras_3[13],tras_3[14],tras_3[15]};


tras_4[]=Translate {lx*2*0, ly*2*4, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_4") ={tras_4[4],tras_4[5],tras_4[6],tras_4[7],tras_4[8],tras_4[9],tras_4[10],tras_4[11]};
Physical Surface("MATRIX_4") ={tras_4[0],tras_4[1],tras_4[2],tras_4[3],tras_4[12],tras_4[13],tras_4[14],tras_4[15]};


tras_5[]=Translate {lx*2*0, ly*2*5, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_5") ={tras_5[4],tras_5[5],tras_5[6],tras_5[7],tras_5[8],tras_5[9],tras_5[10],tras_5[11]};
Physical Surface("MATRIX_5") ={tras_5[0],tras_5[1],tras_5[2],tras_5[3],tras_5[12],tras_5[13],tras_5[14],tras_5[15]};


tras_6[]=Translate {lx*2*0, ly*2*6, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_6") ={tras_6[4],tras_6[5],tras_6[6],tras_6[7],tras_6[8],tras_6[9],tras_6[10],tras_6[11]};
Physical Surface("MATRIX_6") ={tras_6[0],tras_6[1],tras_6[2],tras_6[3],tras_6[12],tras_6[13],tras_6[14],tras_6[15]};


tras_7[]=Translate {lx*2*0, ly*2*7, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_7") ={tras_7[4],tras_7[5],tras_7[6],tras_7[7],tras_7[8],tras_7[9],tras_7[10],tras_7[11]};
Physical Surface("MATRIX_7") ={tras_7[0],tras_7[1],tras_7[2],tras_7[3],tras_7[12],tras_7[13],tras_7[14],tras_7[15]};


tras_8[]=Translate {lx*2*0, ly*2*8, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_8") ={tras_8[4],tras_8[5],tras_8[6],tras_8[7],tras_8[8],tras_8[9],tras_8[10],tras_8[11]};
Physical Surface("MATRIX_8") ={tras_8[0],tras_8[1],tras_8[2],tras_8[3],tras_8[12],tras_8[13],tras_8[14],tras_8[15]};


tras_9[]=Translate {lx*2*0, ly*2*9, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_9") ={tras_9[4],tras_9[5],tras_9[6],tras_9[7],tras_9[8],tras_9[9],tras_9[10],tras_9[11]};
Physical Surface("MATRIX_9") ={tras_9[0],tras_9[1],tras_9[2],tras_9[3],tras_9[12],tras_9[13],tras_9[14],tras_9[15]};


tras_10[]=Translate {lx*2*1, ly*2*0, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_10") ={tras_10[4],tras_10[5],tras_10[6],tras_10[7],tras_10[8],tras_10[9],tras_10[10],tras_10[11]};
Physical Surface("MATRIX_10") ={tras_10[0],tras_10[1],tras_10[2],tras_10[3],tras_10[12],tras_10[13],tras_10[14],tras_10[15]};


tras_11[]=Translate {lx*2*1, ly*2*1, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_11") ={tras_11[4],tras_11[5],tras_11[6],tras_11[7],tras_11[8],tras_11[9],tras_11[10],tras_11[11]};
Physical Surface("MATRIX_11") ={tras_11[0],tras_11[1],tras_11[2],tras_11[3],tras_11[12],tras_11[13],tras_11[14],tras_11[15]};


tras_12[]=Translate {lx*2*1, ly*2*2, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_12") ={tras_12[4],tras_12[5],tras_12[6],tras_12[7],tras_12[8],tras_12[9],tras_12[10],tras_12[11]};
Physical Surface("MATRIX_12") ={tras_12[0],tras_12[1],tras_12[2],tras_12[3],tras_12[12],tras_12[13],tras_12[14],tras_12[15]};


tras_13[]=Translate {lx*2*1, ly*2*3, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_13") ={tras_13[4],tras_13[5],tras_13[6],tras_13[7],tras_13[8],tras_13[9],tras_13[10],tras_13[11]};
Physical Surface("MATRIX_13") ={tras_13[0],tras_13[1],tras_13[2],tras_13[3],tras_13[12],tras_13[13],tras_13[14],tras_13[15]};


tras_14[]=Translate {lx*2*1, ly*2*4, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_14") ={tras_14[4],tras_14[5],tras_14[6],tras_14[7],tras_14[8],tras_14[9],tras_14[10],tras_14[11]};
Physical Surface("MATRIX_14") ={tras_14[0],tras_14[1],tras_14[2],tras_14[3],tras_14[12],tras_14[13],tras_14[14],tras_14[15]};


tras_15[]=Translate {lx*2*1, ly*2*5, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_15") ={tras_15[4],tras_15[5],tras_15[6],tras_15[7],tras_15[8],tras_15[9],tras_15[10],tras_15[11]};
Physical Surface("MATRIX_15") ={tras_15[0],tras_15[1],tras_15[2],tras_15[3],tras_15[12],tras_15[13],tras_15[14],tras_15[15]};


tras_16[]=Translate {lx*2*1, ly*2*6, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_16") ={tras_16[4],tras_16[5],tras_16[6],tras_16[7],tras_16[8],tras_16[9],tras_16[10],tras_16[11]};
Physical Surface("MATRIX_16") ={tras_16[0],tras_16[1],tras_16[2],tras_16[3],tras_16[12],tras_16[13],tras_16[14],tras_16[15]};


tras_17[]=Translate {lx*2*1, ly*2*7, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_17") ={tras_17[4],tras_17[5],tras_17[6],tras_17[7],tras_17[8],tras_17[9],tras_17[10],tras_17[11]};
Physical Surface("MATRIX_17") ={tras_17[0],tras_17[1],tras_17[2],tras_17[3],tras_17[12],tras_17[13],tras_17[14],tras_17[15]};


tras_18[]=Translate {lx*2*1, ly*2*8, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_18") ={tras_18[4],tras_18[5],tras_18[6],tras_18[7],tras_18[8],tras_18[9],tras_18[10],tras_18[11]};
Physical Surface("MATRIX_18") ={tras_18[0],tras_18[1],tras_18[2],tras_18[3],tras_18[12],tras_18[13],tras_18[14],tras_18[15]};


tras_19[]=Translate {lx*2*1, ly*2*9, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_19") ={tras_19[4],tras_19[5],tras_19[6],tras_19[7],tras_19[8],tras_19[9],tras_19[10],tras_19[11]};
Physical Surface("MATRIX_19") ={tras_19[0],tras_19[1],tras_19[2],tras_19[3],tras_19[12],tras_19[13],tras_19[14],tras_19[15]};


tras_20[]=Translate {lx*2*2, ly*2*0, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_20") ={tras_20[4],tras_20[5],tras_20[6],tras_20[7],tras_20[8],tras_20[9],tras_20[10],tras_20[11]};
Physical Surface("MATRIX_20") ={tras_20[0],tras_20[1],tras_20[2],tras_20[3],tras_20[12],tras_20[13],tras_20[14],tras_20[15]};


tras_21[]=Translate {lx*2*2, ly*2*1, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_21") ={tras_21[4],tras_21[5],tras_21[6],tras_21[7],tras_21[8],tras_21[9],tras_21[10],tras_21[11]};
Physical Surface("MATRIX_21") ={tras_21[0],tras_21[1],tras_21[2],tras_21[3],tras_21[12],tras_21[13],tras_21[14],tras_21[15]};


tras_22[]=Translate {lx*2*2, ly*2*2, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_22") ={tras_22[4],tras_22[5],tras_22[6],tras_22[7],tras_22[8],tras_22[9],tras_22[10],tras_22[11]};
Physical Surface("MATRIX_22") ={tras_22[0],tras_22[1],tras_22[2],tras_22[3],tras_22[12],tras_22[13],tras_22[14],tras_22[15]};


tras_23[]=Translate {lx*2*2, ly*2*3, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_23") ={tras_23[4],tras_23[5],tras_23[6],tras_23[7],tras_23[8],tras_23[9],tras_23[10],tras_23[11]};
Physical Surface("MATRIX_23") ={tras_23[0],tras_23[1],tras_23[2],tras_23[3],tras_23[12],tras_23[13],tras_23[14],tras_23[15]};


tras_24[]=Translate {lx*2*2, ly*2*4, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_24") ={tras_24[4],tras_24[5],tras_24[6],tras_24[7],tras_24[8],tras_24[9],tras_24[10],tras_24[11]};
Physical Surface("MATRIX_24") ={tras_24[0],tras_24[1],tras_24[2],tras_24[3],tras_24[12],tras_24[13],tras_24[14],tras_24[15]};


tras_25[]=Translate {lx*2*2, ly*2*5, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_25") ={tras_25[4],tras_25[5],tras_25[6],tras_25[7],tras_25[8],tras_25[9],tras_25[10],tras_25[11]};
Physical Surface("MATRIX_25") ={tras_25[0],tras_25[1],tras_25[2],tras_25[3],tras_25[12],tras_25[13],tras_25[14],tras_25[15]};


tras_26[]=Translate {lx*2*2, ly*2*6, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_26") ={tras_26[4],tras_26[5],tras_26[6],tras_26[7],tras_26[8],tras_26[9],tras_26[10],tras_26[11]};
Physical Surface("MATRIX_26") ={tras_26[0],tras_26[1],tras_26[2],tras_26[3],tras_26[12],tras_26[13],tras_26[14],tras_26[15]};


tras_27[]=Translate {lx*2*2, ly*2*7, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_27") ={tras_27[4],tras_27[5],tras_27[6],tras_27[7],tras_27[8],tras_27[9],tras_27[10],tras_27[11]};
Physical Surface("MATRIX_27") ={tras_27[0],tras_27[1],tras_27[2],tras_27[3],tras_27[12],tras_27[13],tras_27[14],tras_27[15]};


tras_28[]=Translate {lx*2*2, ly*2*8, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_28") ={tras_28[4],tras_28[5],tras_28[6],tras_28[7],tras_28[8],tras_28[9],tras_28[10],tras_28[11]};
Physical Surface("MATRIX_28") ={tras_28[0],tras_28[1],tras_28[2],tras_28[3],tras_28[12],tras_28[13],tras_28[14],tras_28[15]};


tras_29[]=Translate {lx*2*2, ly*2*9, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_29") ={tras_29[4],tras_29[5],tras_29[6],tras_29[7],tras_29[8],tras_29[9],tras_29[10],tras_29[11]};
Physical Surface("MATRIX_29") ={tras_29[0],tras_29[1],tras_29[2],tras_29[3],tras_29[12],tras_29[13],tras_29[14],tras_29[15]};


tras_30[]=Translate {lx*2*3, ly*2*0, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_30") ={tras_30[4],tras_30[5],tras_30[6],tras_30[7],tras_30[8],tras_30[9],tras_30[10],tras_30[11]};
Physical Surface("MATRIX_30") ={tras_30[0],tras_30[1],tras_30[2],tras_30[3],tras_30[12],tras_30[13],tras_30[14],tras_30[15]};


tras_31[]=Translate {lx*2*3, ly*2*1, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_31") ={tras_31[4],tras_31[5],tras_31[6],tras_31[7],tras_31[8],tras_31[9],tras_31[10],tras_31[11]};
Physical Surface("MATRIX_31") ={tras_31[0],tras_31[1],tras_31[2],tras_31[3],tras_31[12],tras_31[13],tras_31[14],tras_31[15]};


tras_32[]=Translate {lx*2*3, ly*2*2, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_32") ={tras_32[4],tras_32[5],tras_32[6],tras_32[7],tras_32[8],tras_32[9],tras_32[10],tras_32[11]};
Physical Surface("MATRIX_32") ={tras_32[0],tras_32[1],tras_32[2],tras_32[3],tras_32[12],tras_32[13],tras_32[14],tras_32[15]};


tras_33[]=Translate {lx*2*3, ly*2*3, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_33") ={tras_33[4],tras_33[5],tras_33[6],tras_33[7],tras_33[8],tras_33[9],tras_33[10],tras_33[11]};
Physical Surface("MATRIX_33") ={tras_33[0],tras_33[1],tras_33[2],tras_33[3],tras_33[12],tras_33[13],tras_33[14],tras_33[15]};


tras_34[]=Translate {lx*2*3, ly*2*4, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_34") ={tras_34[4],tras_34[5],tras_34[6],tras_34[7],tras_34[8],tras_34[9],tras_34[10],tras_34[11]};
Physical Surface("MATRIX_34") ={tras_34[0],tras_34[1],tras_34[2],tras_34[3],tras_34[12],tras_34[13],tras_34[14],tras_34[15]};


tras_35[]=Translate {lx*2*3, ly*2*5, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_35") ={tras_35[4],tras_35[5],tras_35[6],tras_35[7],tras_35[8],tras_35[9],tras_35[10],tras_35[11]};
Physical Surface("MATRIX_35") ={tras_35[0],tras_35[1],tras_35[2],tras_35[3],tras_35[12],tras_35[13],tras_35[14],tras_35[15]};


tras_36[]=Translate {lx*2*3, ly*2*6, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_36") ={tras_36[4],tras_36[5],tras_36[6],tras_36[7],tras_36[8],tras_36[9],tras_36[10],tras_36[11]};
Physical Surface("MATRIX_36") ={tras_36[0],tras_36[1],tras_36[2],tras_36[3],tras_36[12],tras_36[13],tras_36[14],tras_36[15]};


tras_37[]=Translate {lx*2*3, ly*2*7, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_37") ={tras_37[4],tras_37[5],tras_37[6],tras_37[7],tras_37[8],tras_37[9],tras_37[10],tras_37[11]};
Physical Surface("MATRIX_37") ={tras_37[0],tras_37[1],tras_37[2],tras_37[3],tras_37[12],tras_37[13],tras_37[14],tras_37[15]};


tras_38[]=Translate {lx*2*3, ly*2*8, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_38") ={tras_38[4],tras_38[5],tras_38[6],tras_38[7],tras_38[8],tras_38[9],tras_38[10],tras_38[11]};
Physical Surface("MATRIX_38") ={tras_38[0],tras_38[1],tras_38[2],tras_38[3],tras_38[12],tras_38[13],tras_38[14],tras_38[15]};


tras_39[]=Translate {lx*2*3, ly*2*9, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_39") ={tras_39[4],tras_39[5],tras_39[6],tras_39[7],tras_39[8],tras_39[9],tras_39[10],tras_39[11]};
Physical Surface("MATRIX_39") ={tras_39[0],tras_39[1],tras_39[2],tras_39[3],tras_39[12],tras_39[13],tras_39[14],tras_39[15]};


tras_40[]=Translate {lx*2*4, ly*2*0, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_40") ={tras_40[4],tras_40[5],tras_40[6],tras_40[7],tras_40[8],tras_40[9],tras_40[10],tras_40[11]};
Physical Surface("MATRIX_40") ={tras_40[0],tras_40[1],tras_40[2],tras_40[3],tras_40[12],tras_40[13],tras_40[14],tras_40[15]};


tras_41[]=Translate {lx*2*4, ly*2*1, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_41") ={tras_41[4],tras_41[5],tras_41[6],tras_41[7],tras_41[8],tras_41[9],tras_41[10],tras_41[11]};
Physical Surface("MATRIX_41") ={tras_41[0],tras_41[1],tras_41[2],tras_41[3],tras_41[12],tras_41[13],tras_41[14],tras_41[15]};


tras_42[]=Translate {lx*2*4, ly*2*2, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_42") ={tras_42[4],tras_42[5],tras_42[6],tras_42[7],tras_42[8],tras_42[9],tras_42[10],tras_42[11]};
Physical Surface("MATRIX_42") ={tras_42[0],tras_42[1],tras_42[2],tras_42[3],tras_42[12],tras_42[13],tras_42[14],tras_42[15]};


tras_43[]=Translate {lx*2*4, ly*2*3, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_43") ={tras_43[4],tras_43[5],tras_43[6],tras_43[7],tras_43[8],tras_43[9],tras_43[10],tras_43[11]};
Physical Surface("MATRIX_43") ={tras_43[0],tras_43[1],tras_43[2],tras_43[3],tras_43[12],tras_43[13],tras_43[14],tras_43[15]};


tras_44[]=Translate {lx*2*4, ly*2*4, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_44") ={tras_44[4],tras_44[5],tras_44[6],tras_44[7],tras_44[8],tras_44[9],tras_44[10],tras_44[11]};
Physical Surface("MATRIX_44") ={tras_44[0],tras_44[1],tras_44[2],tras_44[3],tras_44[12],tras_44[13],tras_44[14],tras_44[15]};


tras_45[]=Translate {lx*2*4, ly*2*5, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_45") ={tras_45[4],tras_45[5],tras_45[6],tras_45[7],tras_45[8],tras_45[9],tras_45[10],tras_45[11]};
Physical Surface("MATRIX_45") ={tras_45[0],tras_45[1],tras_45[2],tras_45[3],tras_45[12],tras_45[13],tras_45[14],tras_45[15]};


tras_46[]=Translate {lx*2*4, ly*2*6, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_46") ={tras_46[4],tras_46[5],tras_46[6],tras_46[7],tras_46[8],tras_46[9],tras_46[10],tras_46[11]};
Physical Surface("MATRIX_46") ={tras_46[0],tras_46[1],tras_46[2],tras_46[3],tras_46[12],tras_46[13],tras_46[14],tras_46[15]};


tras_47[]=Translate {lx*2*4, ly*2*7, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_47") ={tras_47[4],tras_47[5],tras_47[6],tras_47[7],tras_47[8],tras_47[9],tras_47[10],tras_47[11]};
Physical Surface("MATRIX_47") ={tras_47[0],tras_47[1],tras_47[2],tras_47[3],tras_47[12],tras_47[13],tras_47[14],tras_47[15]};


tras_48[]=Translate {lx*2*4, ly*2*8, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_48") ={tras_48[4],tras_48[5],tras_48[6],tras_48[7],tras_48[8],tras_48[9],tras_48[10],tras_48[11]};
Physical Surface("MATRIX_48") ={tras_48[0],tras_48[1],tras_48[2],tras_48[3],tras_48[12],tras_48[13],tras_48[14],tras_48[15]};


tras_49[]=Translate {lx*2*4, ly*2*9, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_49") ={tras_49[4],tras_49[5],tras_49[6],tras_49[7],tras_49[8],tras_49[9],tras_49[10],tras_49[11]};
Physical Surface("MATRIX_49") ={tras_49[0],tras_49[1],tras_49[2],tras_49[3],tras_49[12],tras_49[13],tras_49[14],tras_49[15]};


tras_50[]=Translate {lx*2*5, ly*2*0, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_50") ={tras_50[4],tras_50[5],tras_50[6],tras_50[7],tras_50[8],tras_50[9],tras_50[10],tras_50[11]};
Physical Surface("MATRIX_50") ={tras_50[0],tras_50[1],tras_50[2],tras_50[3],tras_50[12],tras_50[13],tras_50[14],tras_50[15]};


tras_51[]=Translate {lx*2*5, ly*2*1, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_51") ={tras_51[4],tras_51[5],tras_51[6],tras_51[7],tras_51[8],tras_51[9],tras_51[10],tras_51[11]};
Physical Surface("MATRIX_51") ={tras_51[0],tras_51[1],tras_51[2],tras_51[3],tras_51[12],tras_51[13],tras_51[14],tras_51[15]};


tras_52[]=Translate {lx*2*5, ly*2*2, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_52") ={tras_52[4],tras_52[5],tras_52[6],tras_52[7],tras_52[8],tras_52[9],tras_52[10],tras_52[11]};
Physical Surface("MATRIX_52") ={tras_52[0],tras_52[1],tras_52[2],tras_52[3],tras_52[12],tras_52[13],tras_52[14],tras_52[15]};


tras_53[]=Translate {lx*2*5, ly*2*3, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_53") ={tras_53[4],tras_53[5],tras_53[6],tras_53[7],tras_53[8],tras_53[9],tras_53[10],tras_53[11]};
Physical Surface("MATRIX_53") ={tras_53[0],tras_53[1],tras_53[2],tras_53[3],tras_53[12],tras_53[13],tras_53[14],tras_53[15]};


tras_54[]=Translate {lx*2*5, ly*2*4, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_54") ={tras_54[4],tras_54[5],tras_54[6],tras_54[7],tras_54[8],tras_54[9],tras_54[10],tras_54[11]};
Physical Surface("MATRIX_54") ={tras_54[0],tras_54[1],tras_54[2],tras_54[3],tras_54[12],tras_54[13],tras_54[14],tras_54[15]};


tras_55[]=Translate {lx*2*5, ly*2*5, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_55") ={tras_55[4],tras_55[5],tras_55[6],tras_55[7],tras_55[8],tras_55[9],tras_55[10],tras_55[11]};
Physical Surface("MATRIX_55") ={tras_55[0],tras_55[1],tras_55[2],tras_55[3],tras_55[12],tras_55[13],tras_55[14],tras_55[15]};


tras_56[]=Translate {lx*2*5, ly*2*6, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_56") ={tras_56[4],tras_56[5],tras_56[6],tras_56[7],tras_56[8],tras_56[9],tras_56[10],tras_56[11]};
Physical Surface("MATRIX_56") ={tras_56[0],tras_56[1],tras_56[2],tras_56[3],tras_56[12],tras_56[13],tras_56[14],tras_56[15]};


tras_57[]=Translate {lx*2*5, ly*2*7, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_57") ={tras_57[4],tras_57[5],tras_57[6],tras_57[7],tras_57[8],tras_57[9],tras_57[10],tras_57[11]};
Physical Surface("MATRIX_57") ={tras_57[0],tras_57[1],tras_57[2],tras_57[3],tras_57[12],tras_57[13],tras_57[14],tras_57[15]};


tras_58[]=Translate {lx*2*5, ly*2*8, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_58") ={tras_58[4],tras_58[5],tras_58[6],tras_58[7],tras_58[8],tras_58[9],tras_58[10],tras_58[11]};
Physical Surface("MATRIX_58") ={tras_58[0],tras_58[1],tras_58[2],tras_58[3],tras_58[12],tras_58[13],tras_58[14],tras_58[15]};


tras_59[]=Translate {lx*2*5, ly*2*9, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_59") ={tras_59[4],tras_59[5],tras_59[6],tras_59[7],tras_59[8],tras_59[9],tras_59[10],tras_59[11]};
Physical Surface("MATRIX_59") ={tras_59[0],tras_59[1],tras_59[2],tras_59[3],tras_59[12],tras_59[13],tras_59[14],tras_59[15]};


tras_60[]=Translate {lx*2*6, ly*2*0, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_60") ={tras_60[4],tras_60[5],tras_60[6],tras_60[7],tras_60[8],tras_60[9],tras_60[10],tras_60[11]};
Physical Surface("MATRIX_60") ={tras_60[0],tras_60[1],tras_60[2],tras_60[3],tras_60[12],tras_60[13],tras_60[14],tras_60[15]};


tras_61[]=Translate {lx*2*6, ly*2*1, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_61") ={tras_61[4],tras_61[5],tras_61[6],tras_61[7],tras_61[8],tras_61[9],tras_61[10],tras_61[11]};
Physical Surface("MATRIX_61") ={tras_61[0],tras_61[1],tras_61[2],tras_61[3],tras_61[12],tras_61[13],tras_61[14],tras_61[15]};


tras_62[]=Translate {lx*2*6, ly*2*2, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_62") ={tras_62[4],tras_62[5],tras_62[6],tras_62[7],tras_62[8],tras_62[9],tras_62[10],tras_62[11]};
Physical Surface("MATRIX_62") ={tras_62[0],tras_62[1],tras_62[2],tras_62[3],tras_62[12],tras_62[13],tras_62[14],tras_62[15]};


tras_63[]=Translate {lx*2*6, ly*2*3, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_63") ={tras_63[4],tras_63[5],tras_63[6],tras_63[7],tras_63[8],tras_63[9],tras_63[10],tras_63[11]};
Physical Surface("MATRIX_63") ={tras_63[0],tras_63[1],tras_63[2],tras_63[3],tras_63[12],tras_63[13],tras_63[14],tras_63[15]};


tras_64[]=Translate {lx*2*6, ly*2*4, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_64") ={tras_64[4],tras_64[5],tras_64[6],tras_64[7],tras_64[8],tras_64[9],tras_64[10],tras_64[11]};
Physical Surface("MATRIX_64") ={tras_64[0],tras_64[1],tras_64[2],tras_64[3],tras_64[12],tras_64[13],tras_64[14],tras_64[15]};


tras_65[]=Translate {lx*2*6, ly*2*5, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_65") ={tras_65[4],tras_65[5],tras_65[6],tras_65[7],tras_65[8],tras_65[9],tras_65[10],tras_65[11]};
Physical Surface("MATRIX_65") ={tras_65[0],tras_65[1],tras_65[2],tras_65[3],tras_65[12],tras_65[13],tras_65[14],tras_65[15]};


tras_66[]=Translate {lx*2*6, ly*2*6, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_66") ={tras_66[4],tras_66[5],tras_66[6],tras_66[7],tras_66[8],tras_66[9],tras_66[10],tras_66[11]};
Physical Surface("MATRIX_66") ={tras_66[0],tras_66[1],tras_66[2],tras_66[3],tras_66[12],tras_66[13],tras_66[14],tras_66[15]};


tras_67[]=Translate {lx*2*6, ly*2*7, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_67") ={tras_67[4],tras_67[5],tras_67[6],tras_67[7],tras_67[8],tras_67[9],tras_67[10],tras_67[11]};
Physical Surface("MATRIX_67") ={tras_67[0],tras_67[1],tras_67[2],tras_67[3],tras_67[12],tras_67[13],tras_67[14],tras_67[15]};


tras_68[]=Translate {lx*2*6, ly*2*8, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_68") ={tras_68[4],tras_68[5],tras_68[6],tras_68[7],tras_68[8],tras_68[9],tras_68[10],tras_68[11]};
Physical Surface("MATRIX_68") ={tras_68[0],tras_68[1],tras_68[2],tras_68[3],tras_68[12],tras_68[13],tras_68[14],tras_68[15]};


tras_69[]=Translate {lx*2*6, ly*2*9, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_69") ={tras_69[4],tras_69[5],tras_69[6],tras_69[7],tras_69[8],tras_69[9],tras_69[10],tras_69[11]};
Physical Surface("MATRIX_69") ={tras_69[0],tras_69[1],tras_69[2],tras_69[3],tras_69[12],tras_69[13],tras_69[14],tras_69[15]};


tras_70[]=Translate {lx*2*7, ly*2*0, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_70") ={tras_70[4],tras_70[5],tras_70[6],tras_70[7],tras_70[8],tras_70[9],tras_70[10],tras_70[11]};
Physical Surface("MATRIX_70") ={tras_70[0],tras_70[1],tras_70[2],tras_70[3],tras_70[12],tras_70[13],tras_70[14],tras_70[15]};


tras_71[]=Translate {lx*2*7, ly*2*1, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_71") ={tras_71[4],tras_71[5],tras_71[6],tras_71[7],tras_71[8],tras_71[9],tras_71[10],tras_71[11]};
Physical Surface("MATRIX_71") ={tras_71[0],tras_71[1],tras_71[2],tras_71[3],tras_71[12],tras_71[13],tras_71[14],tras_71[15]};


tras_72[]=Translate {lx*2*7, ly*2*2, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_72") ={tras_72[4],tras_72[5],tras_72[6],tras_72[7],tras_72[8],tras_72[9],tras_72[10],tras_72[11]};
Physical Surface("MATRIX_72") ={tras_72[0],tras_72[1],tras_72[2],tras_72[3],tras_72[12],tras_72[13],tras_72[14],tras_72[15]};


tras_73[]=Translate {lx*2*7, ly*2*3, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_73") ={tras_73[4],tras_73[5],tras_73[6],tras_73[7],tras_73[8],tras_73[9],tras_73[10],tras_73[11]};
Physical Surface("MATRIX_73") ={tras_73[0],tras_73[1],tras_73[2],tras_73[3],tras_73[12],tras_73[13],tras_73[14],tras_73[15]};


tras_74[]=Translate {lx*2*7, ly*2*4, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_74") ={tras_74[4],tras_74[5],tras_74[6],tras_74[7],tras_74[8],tras_74[9],tras_74[10],tras_74[11]};
Physical Surface("MATRIX_74") ={tras_74[0],tras_74[1],tras_74[2],tras_74[3],tras_74[12],tras_74[13],tras_74[14],tras_74[15]};


tras_75[]=Translate {lx*2*7, ly*2*5, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_75") ={tras_75[4],tras_75[5],tras_75[6],tras_75[7],tras_75[8],tras_75[9],tras_75[10],tras_75[11]};
Physical Surface("MATRIX_75") ={tras_75[0],tras_75[1],tras_75[2],tras_75[3],tras_75[12],tras_75[13],tras_75[14],tras_75[15]};


tras_76[]=Translate {lx*2*7, ly*2*6, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_76") ={tras_76[4],tras_76[5],tras_76[6],tras_76[7],tras_76[8],tras_76[9],tras_76[10],tras_76[11]};
Physical Surface("MATRIX_76") ={tras_76[0],tras_76[1],tras_76[2],tras_76[3],tras_76[12],tras_76[13],tras_76[14],tras_76[15]};


tras_77[]=Translate {lx*2*7, ly*2*7, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_77") ={tras_77[4],tras_77[5],tras_77[6],tras_77[7],tras_77[8],tras_77[9],tras_77[10],tras_77[11]};
Physical Surface("MATRIX_77") ={tras_77[0],tras_77[1],tras_77[2],tras_77[3],tras_77[12],tras_77[13],tras_77[14],tras_77[15]};


tras_78[]=Translate {lx*2*7, ly*2*8, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_78") ={tras_78[4],tras_78[5],tras_78[6],tras_78[7],tras_78[8],tras_78[9],tras_78[10],tras_78[11]};
Physical Surface("MATRIX_78") ={tras_78[0],tras_78[1],tras_78[2],tras_78[3],tras_78[12],tras_78[13],tras_78[14],tras_78[15]};


tras_79[]=Translate {lx*2*7, ly*2*9, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_79") ={tras_79[4],tras_79[5],tras_79[6],tras_79[7],tras_79[8],tras_79[9],tras_79[10],tras_79[11]};
Physical Surface("MATRIX_79") ={tras_79[0],tras_79[1],tras_79[2],tras_79[3],tras_79[12],tras_79[13],tras_79[14],tras_79[15]};


tras_80[]=Translate {lx*2*8, ly*2*0, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_80") ={tras_80[4],tras_80[5],tras_80[6],tras_80[7],tras_80[8],tras_80[9],tras_80[10],tras_80[11]};
Physical Surface("MATRIX_80") ={tras_80[0],tras_80[1],tras_80[2],tras_80[3],tras_80[12],tras_80[13],tras_80[14],tras_80[15]};


tras_81[]=Translate {lx*2*8, ly*2*1, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_81") ={tras_81[4],tras_81[5],tras_81[6],tras_81[7],tras_81[8],tras_81[9],tras_81[10],tras_81[11]};
Physical Surface("MATRIX_81") ={tras_81[0],tras_81[1],tras_81[2],tras_81[3],tras_81[12],tras_81[13],tras_81[14],tras_81[15]};


tras_82[]=Translate {lx*2*8, ly*2*2, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_82") ={tras_82[4],tras_82[5],tras_82[6],tras_82[7],tras_82[8],tras_82[9],tras_82[10],tras_82[11]};
Physical Surface("MATRIX_82") ={tras_82[0],tras_82[1],tras_82[2],tras_82[3],tras_82[12],tras_82[13],tras_82[14],tras_82[15]};


tras_83[]=Translate {lx*2*8, ly*2*3, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_83") ={tras_83[4],tras_83[5],tras_83[6],tras_83[7],tras_83[8],tras_83[9],tras_83[10],tras_83[11]};
Physical Surface("MATRIX_83") ={tras_83[0],tras_83[1],tras_83[2],tras_83[3],tras_83[12],tras_83[13],tras_83[14],tras_83[15]};


tras_84[]=Translate {lx*2*8, ly*2*4, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_84") ={tras_84[4],tras_84[5],tras_84[6],tras_84[7],tras_84[8],tras_84[9],tras_84[10],tras_84[11]};
Physical Surface("MATRIX_84") ={tras_84[0],tras_84[1],tras_84[2],tras_84[3],tras_84[12],tras_84[13],tras_84[14],tras_84[15]};


tras_85[]=Translate {lx*2*8, ly*2*5, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_85") ={tras_85[4],tras_85[5],tras_85[6],tras_85[7],tras_85[8],tras_85[9],tras_85[10],tras_85[11]};
Physical Surface("MATRIX_85") ={tras_85[0],tras_85[1],tras_85[2],tras_85[3],tras_85[12],tras_85[13],tras_85[14],tras_85[15]};


tras_86[]=Translate {lx*2*8, ly*2*6, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_86") ={tras_86[4],tras_86[5],tras_86[6],tras_86[7],tras_86[8],tras_86[9],tras_86[10],tras_86[11]};
Physical Surface("MATRIX_86") ={tras_86[0],tras_86[1],tras_86[2],tras_86[3],tras_86[12],tras_86[13],tras_86[14],tras_86[15]};


tras_87[]=Translate {lx*2*8, ly*2*7, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_87") ={tras_87[4],tras_87[5],tras_87[6],tras_87[7],tras_87[8],tras_87[9],tras_87[10],tras_87[11]};
Physical Surface("MATRIX_87") ={tras_87[0],tras_87[1],tras_87[2],tras_87[3],tras_87[12],tras_87[13],tras_87[14],tras_87[15]};


tras_88[]=Translate {lx*2*8, ly*2*8, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_88") ={tras_88[4],tras_88[5],tras_88[6],tras_88[7],tras_88[8],tras_88[9],tras_88[10],tras_88[11]};
Physical Surface("MATRIX_88") ={tras_88[0],tras_88[1],tras_88[2],tras_88[3],tras_88[12],tras_88[13],tras_88[14],tras_88[15]};


tras_89[]=Translate {lx*2*8, ly*2*9, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_89") ={tras_89[4],tras_89[5],tras_89[6],tras_89[7],tras_89[8],tras_89[9],tras_89[10],tras_89[11]};
Physical Surface("MATRIX_89") ={tras_89[0],tras_89[1],tras_89[2],tras_89[3],tras_89[12],tras_89[13],tras_89[14],tras_89[15]};


tras_90[]=Translate {lx*2*9, ly*2*0, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_90") ={tras_90[4],tras_90[5],tras_90[6],tras_90[7],tras_90[8],tras_90[9],tras_90[10],tras_90[11]};
Physical Surface("MATRIX_90") ={tras_90[0],tras_90[1],tras_90[2],tras_90[3],tras_90[12],tras_90[13],tras_90[14],tras_90[15]};


tras_91[]=Translate {lx*2*9, ly*2*1, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_91") ={tras_91[4],tras_91[5],tras_91[6],tras_91[7],tras_91[8],tras_91[9],tras_91[10],tras_91[11]};
Physical Surface("MATRIX_91") ={tras_91[0],tras_91[1],tras_91[2],tras_91[3],tras_91[12],tras_91[13],tras_91[14],tras_91[15]};


tras_92[]=Translate {lx*2*9, ly*2*2, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_92") ={tras_92[4],tras_92[5],tras_92[6],tras_92[7],tras_92[8],tras_92[9],tras_92[10],tras_92[11]};
Physical Surface("MATRIX_92") ={tras_92[0],tras_92[1],tras_92[2],tras_92[3],tras_92[12],tras_92[13],tras_92[14],tras_92[15]};


tras_93[]=Translate {lx*2*9, ly*2*3, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_93") ={tras_93[4],tras_93[5],tras_93[6],tras_93[7],tras_93[8],tras_93[9],tras_93[10],tras_93[11]};
Physical Surface("MATRIX_93") ={tras_93[0],tras_93[1],tras_93[2],tras_93[3],tras_93[12],tras_93[13],tras_93[14],tras_93[15]};


tras_94[]=Translate {lx*2*9, ly*2*4, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_94") ={tras_94[4],tras_94[5],tras_94[6],tras_94[7],tras_94[8],tras_94[9],tras_94[10],tras_94[11]};
Physical Surface("MATRIX_94") ={tras_94[0],tras_94[1],tras_94[2],tras_94[3],tras_94[12],tras_94[13],tras_94[14],tras_94[15]};


tras_95[]=Translate {lx*2*9, ly*2*5, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_95") ={tras_95[4],tras_95[5],tras_95[6],tras_95[7],tras_95[8],tras_95[9],tras_95[10],tras_95[11]};
Physical Surface("MATRIX_95") ={tras_95[0],tras_95[1],tras_95[2],tras_95[3],tras_95[12],tras_95[13],tras_95[14],tras_95[15]};


tras_96[]=Translate {lx*2*9, ly*2*6, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_96") ={tras_96[4],tras_96[5],tras_96[6],tras_96[7],tras_96[8],tras_96[9],tras_96[10],tras_96[11]};
Physical Surface("MATRIX_96") ={tras_96[0],tras_96[1],tras_96[2],tras_96[3],tras_96[12],tras_96[13],tras_96[14],tras_96[15]};


tras_97[]=Translate {lx*2*9, ly*2*7, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_97") ={tras_97[4],tras_97[5],tras_97[6],tras_97[7],tras_97[8],tras_97[9],tras_97[10],tras_97[11]};
Physical Surface("MATRIX_97") ={tras_97[0],tras_97[1],tras_97[2],tras_97[3],tras_97[12],tras_97[13],tras_97[14],tras_97[15]};


tras_98[]=Translate {lx*2*9, ly*2*8, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_98") ={tras_98[4],tras_98[5],tras_98[6],tras_98[7],tras_98[8],tras_98[9],tras_98[10],tras_98[11]};
Physical Surface("MATRIX_98") ={tras_98[0],tras_98[1],tras_98[2],tras_98[3],tras_98[12],tras_98[13],tras_98[14],tras_98[15]};


tras_99[]=Translate {lx*2*9, ly*2*9, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface("FIBER_99") ={tras_99[4],tras_99[5],tras_99[6],tras_99[7],tras_99[8],tras_99[9],tras_99[10],tras_99[11]};
Physical Surface("MATRIX_99") ={tras_99[0],tras_99[1],tras_99[2],tras_99[3],tras_99[12],tras_99[13],tras_99[14],tras_99[15]};

Physical Line(7196) = {3, 39, 777, 782, 1497, 1502, 2217, 2222, 2937, 2942, 3657, 3662, 4377, 4382, 5097, 5102, 5817, 5822, 6537, 6542};
Physical Line(7197) = {653, 658, 1373, 1378, 2093, 2098, 2813, 2818, 3533, 3538, 4253, 4258, 4973, 4978, 5693, 5698, 6413, 6418, 7133, 7138};
Physical Line(7198) = {15, 33, 124, 72, 196, 144, 268, 216, 340, 288, 412, 360, 484, 432, 556, 504, 628, 576, 700, 648};
Physical Line("X0") = {6547, 6495, 6619, 6567, 6691, 6763, 6639, 6711, 6835, 6783, 6907, 6855, 6979, 6927, 7051, 6999, 7123, 7071, 7195, 7143};
Physical Line(7200) = {648, 700, 576, 628, 504, 556, 432, 484, 360, 412, 340, 288, 216, 268, 144, 196, 72, 124, 33, 15};
Physical Line("X1") = {7143, 7195, 7071, 7123, 6999, 7051, 6927, 6979, 6855, 6907, 6783, 6835, 6711, 6763, 6639, 6691, 6567, 6619, 6495, 6547};
