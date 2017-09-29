//
// Basic 2D figure with hole and a symmetry
//

Geometry.CopyMeshingMethod = 1;

lx=3.0;
ly=3.0;
Nx = 20;
Ny = 5;

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
Transfinite Surface {12,14,16};

Recombine Surface {16, 14, 12};

Physical Surface("FIBER_0")  = {14};
Physical Surface("MATRIX_0") = {12, 16};


          tras_1[]=Translate {lx*0, ly*1, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_1")  ={tras_1[1]};
          Physical Surface("MATRIX_1") ={tras_1[0],tras_1[2]};
          

          tras_2[]=Translate {lx*0, ly*2, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_2")  ={tras_2[1]};
          Physical Surface("MATRIX_2") ={tras_2[0],tras_2[2]};
          

          tras_3[]=Translate {lx*0, ly*3, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_3")  ={tras_3[1]};
          Physical Surface("MATRIX_3") ={tras_3[0],tras_3[2]};
          

          tras_4[]=Translate {lx*0, ly*4, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_4")  ={tras_4[1]};
          Physical Surface("MATRIX_4") ={tras_4[0],tras_4[2]};
          

          tras_5[]=Translate {lx*0, ly*5, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_5")  ={tras_5[1]};
          Physical Surface("MATRIX_5") ={tras_5[0],tras_5[2]};
          

          tras_6[]=Translate {lx*0, ly*6, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_6")  ={tras_6[1]};
          Physical Surface("MATRIX_6") ={tras_6[0],tras_6[2]};
          

          tras_7[]=Translate {lx*0, ly*7, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_7")  ={tras_7[1]};
          Physical Surface("MATRIX_7") ={tras_7[0],tras_7[2]};
          

          tras_8[]=Translate {lx*0, ly*8, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_8")  ={tras_8[1]};
          Physical Surface("MATRIX_8") ={tras_8[0],tras_8[2]};
          

          tras_9[]=Translate {lx*0, ly*9, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_9")  ={tras_9[1]};
          Physical Surface("MATRIX_9") ={tras_9[0],tras_9[2]};
          

          tras_10[]=Translate {lx*1, ly*0, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_10")  ={tras_10[1]};
          Physical Surface("MATRIX_10") ={tras_10[0],tras_10[2]};
          

          tras_11[]=Translate {lx*1, ly*1, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_11")  ={tras_11[1]};
          Physical Surface("MATRIX_11") ={tras_11[0],tras_11[2]};
          

          tras_12[]=Translate {lx*1, ly*2, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_12")  ={tras_12[1]};
          Physical Surface("MATRIX_12") ={tras_12[0],tras_12[2]};
          

          tras_13[]=Translate {lx*1, ly*3, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_13")  ={tras_13[1]};
          Physical Surface("MATRIX_13") ={tras_13[0],tras_13[2]};
          

          tras_14[]=Translate {lx*1, ly*4, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_14")  ={tras_14[1]};
          Physical Surface("MATRIX_14") ={tras_14[0],tras_14[2]};
          

          tras_15[]=Translate {lx*1, ly*5, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_15")  ={tras_15[1]};
          Physical Surface("MATRIX_15") ={tras_15[0],tras_15[2]};
          

          tras_16[]=Translate {lx*1, ly*6, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_16")  ={tras_16[1]};
          Physical Surface("MATRIX_16") ={tras_16[0],tras_16[2]};
          

          tras_17[]=Translate {lx*1, ly*7, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_17")  ={tras_17[1]};
          Physical Surface("MATRIX_17") ={tras_17[0],tras_17[2]};
          

          tras_18[]=Translate {lx*1, ly*8, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_18")  ={tras_18[1]};
          Physical Surface("MATRIX_18") ={tras_18[0],tras_18[2]};
          

          tras_19[]=Translate {lx*1, ly*9, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_19")  ={tras_19[1]};
          Physical Surface("MATRIX_19") ={tras_19[0],tras_19[2]};
          

          tras_20[]=Translate {lx*2, ly*0, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_20")  ={tras_20[1]};
          Physical Surface("MATRIX_20") ={tras_20[0],tras_20[2]};
          

          tras_21[]=Translate {lx*2, ly*1, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_21")  ={tras_21[1]};
          Physical Surface("MATRIX_21") ={tras_21[0],tras_21[2]};
          

          tras_22[]=Translate {lx*2, ly*2, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_22")  ={tras_22[1]};
          Physical Surface("MATRIX_22") ={tras_22[0],tras_22[2]};
          

          tras_23[]=Translate {lx*2, ly*3, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_23")  ={tras_23[1]};
          Physical Surface("MATRIX_23") ={tras_23[0],tras_23[2]};
          

          tras_24[]=Translate {lx*2, ly*4, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_24")  ={tras_24[1]};
          Physical Surface("MATRIX_24") ={tras_24[0],tras_24[2]};
          

          tras_25[]=Translate {lx*2, ly*5, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_25")  ={tras_25[1]};
          Physical Surface("MATRIX_25") ={tras_25[0],tras_25[2]};
          

          tras_26[]=Translate {lx*2, ly*6, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_26")  ={tras_26[1]};
          Physical Surface("MATRIX_26") ={tras_26[0],tras_26[2]};
          

          tras_27[]=Translate {lx*2, ly*7, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_27")  ={tras_27[1]};
          Physical Surface("MATRIX_27") ={tras_27[0],tras_27[2]};
          

          tras_28[]=Translate {lx*2, ly*8, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_28")  ={tras_28[1]};
          Physical Surface("MATRIX_28") ={tras_28[0],tras_28[2]};
          

          tras_29[]=Translate {lx*2, ly*9, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_29")  ={tras_29[1]};
          Physical Surface("MATRIX_29") ={tras_29[0],tras_29[2]};
          

          tras_30[]=Translate {lx*3, ly*0, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_30")  ={tras_30[1]};
          Physical Surface("MATRIX_30") ={tras_30[0],tras_30[2]};
          

          tras_31[]=Translate {lx*3, ly*1, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_31")  ={tras_31[1]};
          Physical Surface("MATRIX_31") ={tras_31[0],tras_31[2]};
          

          tras_32[]=Translate {lx*3, ly*2, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_32")  ={tras_32[1]};
          Physical Surface("MATRIX_32") ={tras_32[0],tras_32[2]};
          

          tras_33[]=Translate {lx*3, ly*3, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_33")  ={tras_33[1]};
          Physical Surface("MATRIX_33") ={tras_33[0],tras_33[2]};
          

          tras_34[]=Translate {lx*3, ly*4, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_34")  ={tras_34[1]};
          Physical Surface("MATRIX_34") ={tras_34[0],tras_34[2]};
          

          tras_35[]=Translate {lx*3, ly*5, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_35")  ={tras_35[1]};
          Physical Surface("MATRIX_35") ={tras_35[0],tras_35[2]};
          

          tras_36[]=Translate {lx*3, ly*6, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_36")  ={tras_36[1]};
          Physical Surface("MATRIX_36") ={tras_36[0],tras_36[2]};
          

          tras_37[]=Translate {lx*3, ly*7, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_37")  ={tras_37[1]};
          Physical Surface("MATRIX_37") ={tras_37[0],tras_37[2]};
          

          tras_38[]=Translate {lx*3, ly*8, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_38")  ={tras_38[1]};
          Physical Surface("MATRIX_38") ={tras_38[0],tras_38[2]};
          

          tras_39[]=Translate {lx*3, ly*9, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_39")  ={tras_39[1]};
          Physical Surface("MATRIX_39") ={tras_39[0],tras_39[2]};
          

          tras_40[]=Translate {lx*4, ly*0, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_40")  ={tras_40[1]};
          Physical Surface("MATRIX_40") ={tras_40[0],tras_40[2]};
          

          tras_41[]=Translate {lx*4, ly*1, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_41")  ={tras_41[1]};
          Physical Surface("MATRIX_41") ={tras_41[0],tras_41[2]};
          

          tras_42[]=Translate {lx*4, ly*2, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_42")  ={tras_42[1]};
          Physical Surface("MATRIX_42") ={tras_42[0],tras_42[2]};
          

          tras_43[]=Translate {lx*4, ly*3, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_43")  ={tras_43[1]};
          Physical Surface("MATRIX_43") ={tras_43[0],tras_43[2]};
          

          tras_44[]=Translate {lx*4, ly*4, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_44")  ={tras_44[1]};
          Physical Surface("MATRIX_44") ={tras_44[0],tras_44[2]};
          

          tras_45[]=Translate {lx*4, ly*5, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_45")  ={tras_45[1]};
          Physical Surface("MATRIX_45") ={tras_45[0],tras_45[2]};
          

          tras_46[]=Translate {lx*4, ly*6, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_46")  ={tras_46[1]};
          Physical Surface("MATRIX_46") ={tras_46[0],tras_46[2]};
          

          tras_47[]=Translate {lx*4, ly*7, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_47")  ={tras_47[1]};
          Physical Surface("MATRIX_47") ={tras_47[0],tras_47[2]};
          

          tras_48[]=Translate {lx*4, ly*8, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_48")  ={tras_48[1]};
          Physical Surface("MATRIX_48") ={tras_48[0],tras_48[2]};
          

          tras_49[]=Translate {lx*4, ly*9, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_49")  ={tras_49[1]};
          Physical Surface("MATRIX_49") ={tras_49[0],tras_49[2]};
          

          tras_50[]=Translate {lx*5, ly*0, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_50")  ={tras_50[1]};
          Physical Surface("MATRIX_50") ={tras_50[0],tras_50[2]};
          

          tras_51[]=Translate {lx*5, ly*1, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_51")  ={tras_51[1]};
          Physical Surface("MATRIX_51") ={tras_51[0],tras_51[2]};
          

          tras_52[]=Translate {lx*5, ly*2, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_52")  ={tras_52[1]};
          Physical Surface("MATRIX_52") ={tras_52[0],tras_52[2]};
          

          tras_53[]=Translate {lx*5, ly*3, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_53")  ={tras_53[1]};
          Physical Surface("MATRIX_53") ={tras_53[0],tras_53[2]};
          

          tras_54[]=Translate {lx*5, ly*4, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_54")  ={tras_54[1]};
          Physical Surface("MATRIX_54") ={tras_54[0],tras_54[2]};
          

          tras_55[]=Translate {lx*5, ly*5, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_55")  ={tras_55[1]};
          Physical Surface("MATRIX_55") ={tras_55[0],tras_55[2]};
          

          tras_56[]=Translate {lx*5, ly*6, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_56")  ={tras_56[1]};
          Physical Surface("MATRIX_56") ={tras_56[0],tras_56[2]};
          

          tras_57[]=Translate {lx*5, ly*7, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_57")  ={tras_57[1]};
          Physical Surface("MATRIX_57") ={tras_57[0],tras_57[2]};
          

          tras_58[]=Translate {lx*5, ly*8, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_58")  ={tras_58[1]};
          Physical Surface("MATRIX_58") ={tras_58[0],tras_58[2]};
          

          tras_59[]=Translate {lx*5, ly*9, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_59")  ={tras_59[1]};
          Physical Surface("MATRIX_59") ={tras_59[0],tras_59[2]};
          

          tras_60[]=Translate {lx*6, ly*0, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_60")  ={tras_60[1]};
          Physical Surface("MATRIX_60") ={tras_60[0],tras_60[2]};
          

          tras_61[]=Translate {lx*6, ly*1, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_61")  ={tras_61[1]};
          Physical Surface("MATRIX_61") ={tras_61[0],tras_61[2]};
          

          tras_62[]=Translate {lx*6, ly*2, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_62")  ={tras_62[1]};
          Physical Surface("MATRIX_62") ={tras_62[0],tras_62[2]};
          

          tras_63[]=Translate {lx*6, ly*3, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_63")  ={tras_63[1]};
          Physical Surface("MATRIX_63") ={tras_63[0],tras_63[2]};
          

          tras_64[]=Translate {lx*6, ly*4, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_64")  ={tras_64[1]};
          Physical Surface("MATRIX_64") ={tras_64[0],tras_64[2]};
          

          tras_65[]=Translate {lx*6, ly*5, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_65")  ={tras_65[1]};
          Physical Surface("MATRIX_65") ={tras_65[0],tras_65[2]};
          

          tras_66[]=Translate {lx*6, ly*6, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_66")  ={tras_66[1]};
          Physical Surface("MATRIX_66") ={tras_66[0],tras_66[2]};
          

          tras_67[]=Translate {lx*6, ly*7, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_67")  ={tras_67[1]};
          Physical Surface("MATRIX_67") ={tras_67[0],tras_67[2]};
          

          tras_68[]=Translate {lx*6, ly*8, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_68")  ={tras_68[1]};
          Physical Surface("MATRIX_68") ={tras_68[0],tras_68[2]};
          

          tras_69[]=Translate {lx*6, ly*9, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_69")  ={tras_69[1]};
          Physical Surface("MATRIX_69") ={tras_69[0],tras_69[2]};
          

          tras_70[]=Translate {lx*7, ly*0, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_70")  ={tras_70[1]};
          Physical Surface("MATRIX_70") ={tras_70[0],tras_70[2]};
          

          tras_71[]=Translate {lx*7, ly*1, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_71")  ={tras_71[1]};
          Physical Surface("MATRIX_71") ={tras_71[0],tras_71[2]};
          

          tras_72[]=Translate {lx*7, ly*2, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_72")  ={tras_72[1]};
          Physical Surface("MATRIX_72") ={tras_72[0],tras_72[2]};
          

          tras_73[]=Translate {lx*7, ly*3, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_73")  ={tras_73[1]};
          Physical Surface("MATRIX_73") ={tras_73[0],tras_73[2]};
          

          tras_74[]=Translate {lx*7, ly*4, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_74")  ={tras_74[1]};
          Physical Surface("MATRIX_74") ={tras_74[0],tras_74[2]};
          

          tras_75[]=Translate {lx*7, ly*5, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_75")  ={tras_75[1]};
          Physical Surface("MATRIX_75") ={tras_75[0],tras_75[2]};
          

          tras_76[]=Translate {lx*7, ly*6, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_76")  ={tras_76[1]};
          Physical Surface("MATRIX_76") ={tras_76[0],tras_76[2]};
          

          tras_77[]=Translate {lx*7, ly*7, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_77")  ={tras_77[1]};
          Physical Surface("MATRIX_77") ={tras_77[0],tras_77[2]};
          

          tras_78[]=Translate {lx*7, ly*8, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_78")  ={tras_78[1]};
          Physical Surface("MATRIX_78") ={tras_78[0],tras_78[2]};
          

          tras_79[]=Translate {lx*7, ly*9, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_79")  ={tras_79[1]};
          Physical Surface("MATRIX_79") ={tras_79[0],tras_79[2]};
          

          tras_80[]=Translate {lx*8, ly*0, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_80")  ={tras_80[1]};
          Physical Surface("MATRIX_80") ={tras_80[0],tras_80[2]};
          

          tras_81[]=Translate {lx*8, ly*1, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_81")  ={tras_81[1]};
          Physical Surface("MATRIX_81") ={tras_81[0],tras_81[2]};
          

          tras_82[]=Translate {lx*8, ly*2, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_82")  ={tras_82[1]};
          Physical Surface("MATRIX_82") ={tras_82[0],tras_82[2]};
          

          tras_83[]=Translate {lx*8, ly*3, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_83")  ={tras_83[1]};
          Physical Surface("MATRIX_83") ={tras_83[0],tras_83[2]};
          

          tras_84[]=Translate {lx*8, ly*4, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_84")  ={tras_84[1]};
          Physical Surface("MATRIX_84") ={tras_84[0],tras_84[2]};
          

          tras_85[]=Translate {lx*8, ly*5, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_85")  ={tras_85[1]};
          Physical Surface("MATRIX_85") ={tras_85[0],tras_85[2]};
          

          tras_86[]=Translate {lx*8, ly*6, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_86")  ={tras_86[1]};
          Physical Surface("MATRIX_86") ={tras_86[0],tras_86[2]};
          

          tras_87[]=Translate {lx*8, ly*7, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_87")  ={tras_87[1]};
          Physical Surface("MATRIX_87") ={tras_87[0],tras_87[2]};
          

          tras_88[]=Translate {lx*8, ly*8, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_88")  ={tras_88[1]};
          Physical Surface("MATRIX_88") ={tras_88[0],tras_88[2]};
          

          tras_89[]=Translate {lx*8, ly*9, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_89")  ={tras_89[1]};
          Physical Surface("MATRIX_89") ={tras_89[0],tras_89[2]};
          

          tras_90[]=Translate {lx*9, ly*0, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_90")  ={tras_90[1]};
          Physical Surface("MATRIX_90") ={tras_90[0],tras_90[2]};
          

          tras_91[]=Translate {lx*9, ly*1, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_91")  ={tras_91[1]};
          Physical Surface("MATRIX_91") ={tras_91[0],tras_91[2]};
          

          tras_92[]=Translate {lx*9, ly*2, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_92")  ={tras_92[1]};
          Physical Surface("MATRIX_92") ={tras_92[0],tras_92[2]};
          

          tras_93[]=Translate {lx*9, ly*3, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_93")  ={tras_93[1]};
          Physical Surface("MATRIX_93") ={tras_93[0],tras_93[2]};
          

          tras_94[]=Translate {lx*9, ly*4, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_94")  ={tras_94[1]};
          Physical Surface("MATRIX_94") ={tras_94[0],tras_94[2]};
          

          tras_95[]=Translate {lx*9, ly*5, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_95")  ={tras_95[1]};
          Physical Surface("MATRIX_95") ={tras_95[0],tras_95[2]};
          

          tras_96[]=Translate {lx*9, ly*6, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_96")  ={tras_96[1]};
          Physical Surface("MATRIX_96") ={tras_96[0],tras_96[2]};
          

          tras_97[]=Translate {lx*9, ly*7, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_97")  ={tras_97[1]};
          Physical Surface("MATRIX_97") ={tras_97[0],tras_97[2]};
          

          tras_98[]=Translate {lx*9, ly*8, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_98")  ={tras_98[1]};
          Physical Surface("MATRIX_98") ={tras_98[0],tras_98[2]};
          

          tras_99[]=Translate {lx*9, ly*9, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface("FIBER_99")  ={tras_99[1]};
          Physical Surface("MATRIX_99") ={tras_99[0],tras_99[2]};
          
Physical Line("Y0") = {4, 165, 296, 427, 558, 689, 820, 951, 1082, 1213};
Physical Line("Y1") = {141, 274, 405, 536, 667, 798, 929, 1060, 1191, 1322};
Physical Line("X0") = {3, 2, 1, 31, 23, 18, 46, 38, 33, 53, 61, 48, 76, 68, 63, 83, 91, 78, 106, 93, 98, 121, 113, 108, 136, 128, 123, 151, 143, 138};
Physical Line("X1") = {1212, 1208, 1203, 1226, 1222, 1217, 1235, 1239, 1230, 1252, 1248, 1243, 1265, 1261, 1278, 1256, 1274, 1269, 1287, 1291, 1304, 1282, 1300, 1295, 1317, 1313, 1308, 1330, 1326, 1321};
