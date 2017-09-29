#!/bin/bash

NX=10
NY=10
file="struct_tr_"$NX"_"$NY".geo"
cp fiber_tr_2d.geo $file

if [ -e "materials.dat" ]; then
   rm materials.dat
fi

for i in $(seq 0 $((NX-1))); do
  for j in $(seq 0 $((NY-1))); do
  
    num=$(($i*NX+$j))

    if [ $num -ne 0 ]; then
    
    echo "
          tras_$num[]=Translate {lx*$i, ly*$j, 0} {
          Duplicata { Surface{12,14,16}; }
          };
          Physical Surface(\"FIBER_$num\")  ={tras_$num[1]};
          Physical Surface(\"MATRIX_$num\") ={tras_$num[0],tras_$num[2]};
          " >> $file
    fi

    echo "MATRIX_$num TYPE00 E=1.0e6 v=0.3"  >> materials.dat
    echo "FIBER_$num  TYPE00 E=1.0e7 v=0.3"  >> materials.dat
  
  done
done
