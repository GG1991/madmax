#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 
EXECUTABLE="/home/guido/codes/sputnik/macro/macro"

if [ -d "direct/force_x" ]; then
  rm -rf direct/force_x/*
else
  mkdir direct/force_x
fi

$MPIEXEC -np 1 $EXECUTABLE  \
    -boundary "X0 11 0 0","X1 11 1 0" \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3"\
    -function "0 2 0.0 0.0 1.0 0.0","1 2 0.0 0.0 1.0 0.01" \
    -mesh direct/direct_10_10_replaced.msh \
    -dim 2 \
    -normal \
    -tf 1.0 \
    -dt 0.2 \
    -pc_type lu \
    -nnz_factor 3 \
    -print_pvtu

mv macro_* direct/force_x/.

if [ -d "direct/force_y" ]; then
  rm -rf direct/force_y/*
else
  mkdir direct/force_y
fi

$MPIEXEC -np 1 $EXECUTABLE  \
    -boundary "X0 11 0 0","X1 11 0 1" \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3"\
    -function "0 2 0.0 0.0 1.0 0.0","1 2 0.0 0.0 1.0 -0.01" \
    -mesh direct/direct_10_10_replaced.msh \
    -dim 2 \
    -normal 1 \
    -tf 1.0 \
    -dt 0.2 \
    -pc_type lu \
    -print_pvtu

mv macro_* direct/force_y/.
