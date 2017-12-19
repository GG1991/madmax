#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 

#$MPIEXEC -np 1 xterm -e gdb -x file_macro.gdb --args ../../macro/macro \
$MPIEXEC -np 1 ../../macro/macro \
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
    -print_pvtu \
: -np 1 ../../micro/micro \
    -coupl \
    -struct_n 75,75 \
    -dim 2 \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3" \
    -micro_struct "fiber_cilin 3.0 3.0 1 1 0.75 0.0 0.0" \
    -pc_type jacobi \
    -ksp_type cg \
    -homo_us \
    -print_pvtu

#xterm -e gdb -x file_micro.gdb --args
#: -np 1 xterm -e gdb -x file_micro.gdb --args ../../micro/micro \
#: -np 1 xterm -e gdb -x file_micro.gdb --args ../../micro/micro \
