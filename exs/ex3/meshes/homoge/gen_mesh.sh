#!/bin/bash

for i in $(seq 2 20); do

  echo "doing mesh $i x $i"
  m4 -Dn_m4=$i homog.geo.m4 > homog.geo
  gmsh -2 homog.geo > /tmp/null
  if [ -e homog.msh ]; then
     mv homog.msh homog_$i.msh
     echo "ok"
  else
     echo "fail"
  fi

done

if [ -e homog.geo ]; then
   rm  homog.geo
fi
