#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 

if [ $# -eq 0 ]; then
  NP=1;
else
  NP=$1;
fi

#$MPIEXEC -np $NP ../../macro/macro \
$MPIEXEC -np $NP xterm -e gdb -x file_macro.gdb --args ../../macro/macro \
    -coupl \
    -boundary "X0 11 0 0","X1 11 1 0" \
    -material "MICRO MAT_MICRO" \
    -function "0 2 0.0 0.0 1.0 0.0","1 2 0.0 0.0 1.0 0.001" \
    -mesh cube_2d.msh \
    -dim 2 \
    -normal \
    -tf 1.0 \
    -dt 0.2 \
    -ksp_type cg \
    -pc_type jacobi \
    -part_geom \
    -nl_max_its 2 \
    -eps_nev 2 \
    -options_left 0   \
: -np $NP xterm -e gdb -x file_micro.gdb --args ../../micro/micro \
    -coupl \
    -struct_n 75,75 \
    -dim 2 \
    -material "matrix mat_elastic 1.0e7 1.0e6 0.3","fiber mat_elastic 1.0e7 1.0e7 0.3" \
    -micro_struct "fiber_cilin 3.0 3.0 1 1 0.75 0.0 0.0" \
    -pc_type jacobi \
    -ksp_type cg \
    -homo_us \
    -options_left  0
