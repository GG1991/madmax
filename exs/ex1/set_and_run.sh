#!/bin/bash
if [ $# -ne 2 ];then
  echo "2 arguments should be given <N> problem size <NP> number of processes"
  exit 1
fi
N=$1
NP=$2

cd ../../meshes/cube_unif
m4 -DN_M4=$N cube.geo.m4 > cube.geo
printf "generating mesh "
gmsh -3 cube.geo 2&> null
cd ../../exs/ex1
printf "0K\n"

./coupling.sh -log_trace
