#!/bin/bash

NX=2
NY=2
L=3.0

file="rve_"$NX"_"$NY".geo"
cp rve_5.geo $file

for i in $(seq 0 $((NY-1))); do
  for j in $(seq 0 $((NX-1))); do
        
    num=$(($i*NX+$j))

    if [ $num -ne 0 ]; then

        lx=$(echo "$j*$L"| bc)
        ly=$(echo "$i*$L"| bc)

	echo "
        trans_$num[] = Translate {$lx, $ly, 0} {
          Duplicata { Surface{10, 12}; }
        };
        
        Physical Surface(\"FIBER_$num\")  = {trans_$num[0]};
        Physical Surface(\"MATRIX_$num\") = {trans_$num[1]};
        " >> $file
    fi

  done
done

