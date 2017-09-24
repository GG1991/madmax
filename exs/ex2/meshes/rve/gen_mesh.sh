#!/bin/bash

for i in $(seq 2 15); do

  echo "doing mesh $i x $i"
  m4 -Dn_m4=$i rve.geo.m4 > rve.geo
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
