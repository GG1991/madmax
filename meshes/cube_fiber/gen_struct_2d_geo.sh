#!/bin/bash

NX=10
NY=10
file="struct_fiber_2d_"$NX"_"$NY".geo"
cp struct_fiber_2d.geo $file

for i in $(seq 0 $((NX-1))); do
  for j in $(seq 0 $((NY-1))); do
  
    num=$(($i*NX+$j))

    if [ $num -ne 0 ]; then
    
    echo "
tras_$num[]=Translate {lx*2*$i, ly*2*$j, 0} {
Duplicata { Surface{29, 20, 53, 62, 67, 58, 34, 25, 10, 16, 40, 49, 11, 8, 35, 44}; }
};
Physical Surface(\"FIBER_$num\") ={tras_$num[4],tras_$num[5],tras_$num[6],tras_$num[7],tras_$num[8],tras_$num[9],tras_$num[10],tras_$num[11]};
Physical Surface(\"MATRIX_$num\") ={tras_$num[0],tras_$num[1],tras_$num[2],tras_$num[3],tras_$num[12],tras_$num[13],tras_$num[14],tras_$num[15]};
" >> $file

    fi
  
  done
done
