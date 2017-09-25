#!/bin/bash

N=11

for i in $(seq 1 15); do

  n=$(echo "$i*$N"| bc)
  echo "doing mesh $n x $n size $l x $l"

  m4 -Dn_m4=$n rve.geo.m4 > rve.geo
  gmsh -2 rve.geo > /tmp/null
  if [ -e rve.msh ]; then
     mv rve.msh rve_$i.msh
     echo "ok"
  else
     echo "fail"
  fi

done

if [ -e rve.geo ]; then
   rm  rve.geo
fi
