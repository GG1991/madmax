#!/bin/bash
if [ $# -ne 1 ];then
  echo "One argument should be given <N> problem size" 
  exit 1
fi
N=$1

cd ../../meshes/cube_unif
m4 -DN_M4=$N cube.geo.m4 > cube.geo
printf "generating mesh "
gmsh -3 cube.geo 2&> null
cd ../../exs/ex1
printf "0K\n"

