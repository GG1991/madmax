#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 
MACRO="/home/guido/codes/sputnik/macro/macro"
MICRO="/home/guido/codes/sputnik/micro/micro"

if [ -d "homog/ts/force_x" ]; then
  rm -rf homog/ts/force_x/*
else
  mkdir homog/ts/force_x
fi

$MPIEXEC -np 1 $MACRO \
    -coupl \
    -boundary "X0 11 0 0","X1 11 1 0" \
    -material "MICRO MAT_MICRO" \
    -function "0 2 0.0 0.0 1.0 0.0","1 2 0.0 0.0 1.0 0.01" \
    -mesh homog/cube_2d.msh \
    -dim 2 \
    -normal \
    -tf 1.0 \
    -dt 0.2 \
    -pc_type lu \
    -part_geom \
    -nl_max_its 2 \
    -eps_nev 2 \
    -print_pvtu \
: -np 1 $MICRO \
    -coupl \
    -struct_n 75,75 \
    -dim 2 \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3" \
    -micro_struct "fiber_line 3.0 3.0 1 1 0.0 0.4 1.0 0.0" \
    -pc_type lu \
    -homo_ts \
    -print_pvtu 

mv macro_* micro_* homog/ts/force_x/.


if [ -d "homog/ts/force_y" ]; then
  rm -rf homog/ts/force_y/*
else
  mkdir homog/ts/force_y
fi

$MPIEXEC -np 1 $MACRO \
    -coupl \
    -boundary "X0 11 0 0","X1 11 0 1" \
    -material "MICRO MAT_MICRO" \
    -function "0 2 0.0 0.0 1.0 0.0","1 2 0.0 0.0 1.0 -0.01" \
    -mesh homog/cube_2d.msh \
    -dim 2 \
    -normal \
    -tf 1.0 \
    -dt 0.2 \
    -pc_type lu \
    -part_geom \
    -nl_max_its 2 \
    -eps_nev 2 \
    -print_pvtu \
: -np 1 $MICRO \
    -coupl \
    -struct_n 75,75 \
    -dim 2 \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3" \
    -micro_struct "fiber_line 3.0 3.0 1 1 0.0 0.4 1.0 0.0" \
    -pc_type lu \
    -homo_ts \
    -print_pvtu 

mv macro_* micro_* homog/ts/force_y/.
