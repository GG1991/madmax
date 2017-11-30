#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 
EXECUTABLE="../../../macro/macro"

if [ $# -eq 0 ]; then
  NP=1;
else
  NP=$1;
fi

#$MPIEXEC -np $NP xterm -e gdb -x file.gdb --args $EXECUTABLE  \
$MPIEXEC -np $NP $EXECUTABLE  \
    -boundary "X0 11 0 0","X1 11 1 0" \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e6 0.3"\
    -function "0 2 0.0 0.0 1.0 0.0","1 2 0.0 0.0 1.0 0.01" \
    -mesh     direct_10_10_replaced.msh \
    -dim 2            \
    -normal           \
    -tf 1.0           \
    -dt 0.2           \
    -pc_type jacobi   \
    -nnz_factor 3     \
    -ksp_type cg      \
    -ksp_atol 1.0e-24 \
    -ksp_dtol 1.0e-10 \
    -ksp_rtol 1.0e-17 \
    -eps_nev  2       \
    -print_vtu        \
    -part_geom        \
    -options_left 0

#-pc_type  jacobi lu \
#-ksp_type cg        \
#-ksp_atol 1.0e-22   \
#-boundary_2 "X1 001 0 0 1"
#-ksp_type cg        \
#-print_petsc        \
