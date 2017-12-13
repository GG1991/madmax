#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 
EXECUTABLE="../../../macro/macro"

if [ $# -eq 0 ]; then
  NP=1;
else
  NP=$1;
fi

if [ -d "force_x" ]; then
  rm -rf force_x/*
else
  mkdir force_x
fi

$MPIEXEC -np $NP $EXECUTABLE  \
    -boundary "X0 11 0 0","X1 11 1 0" \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e6 0.3"\
    -function "0 2 0.0 0.0 1.0 0.0","1 2 0.0 0.0 1.0 0.01" \
    -mesh direct_10_10_replaced.msh \
    -dim 2 \
    -normal \
    -tf 1.0 \
    -dt 0.2 \
    -pc_type jacobi \
    -nnz_factor 3 \
    -ksp_type cg \
    -ksp_atol 1.0e-14 \
    -ksp_rtol 1.0e-13 \
    -print_vtu \
    -options_left 0

mv macro_* force_x/.

if [ -d "force_y" ]; then
  rm -rf force_y/*
else
  mkdir force_y
fi

$MPIEXEC -np $NP $EXECUTABLE  \
    -boundary "X0 11 0 0","X1 11 0 1" \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e6 0.3"\
    -function "0 2 0.0 0.0 1.0 0.0","1 2 0.0 0.0 1.0 -0.01" \
    -mesh direct_10_10_replaced.msh \
    -dim 2 \
    -normal 1 \
    -tf 1.0 \
    -dt 0.2 \
    -pc_type jacobi \
    -ksp_type cg \
    -ksp_atol 1.0e-14 \
    -ksp_rtol 1.0e-13 \
    -print_vtu \
    -options_left 0

mv macro_* force_y/.

#-pc_type  jacobi lu \
#-ksp_type cg        \
#-ksp_atol 1.0e-22   \
#-boundary_2 "X1 001 0 0 1"
#-ksp_type cg        \
#-print_petsc        \
