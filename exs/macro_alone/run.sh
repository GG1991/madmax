#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 

if [ $# -eq 0 ]; then
  NP=1;
else
  NP=$1;
fi

#$MPIEXEC -np $NP xterm -e gdb -x file.gdb --args ../../macro/macro \
$MPIEXEC -np $NP ../../macro/macro \
    -boundary "X0 11 0 0","X1 11 1 0" \
    -material "MATRIX TYPE_0 1.0e7 1.0e6 0.3","MICRO TYPE_1" \
    -function "0 2 0.0 0.0 1.0 0.0","1 2 0.0 0.0 1.0 0.001" \
    -mesh cube_2d.msh \
    -dim 2            \
    -mesh_gmsh        \
    -eigen            \
    -tf 1.0           \
    -dt 1.0           \
    -ksp_type cg      \
    -pc_type jacobi   \
    -part_geom        \
    -nr_max_its  2    \
    -eps_nev     2    \
    -eps_type    power  \
    -print_vtu        \
    -print_petsc      \
    -options_left 0

#xterm -e gdb -x file.gdb --args 
#-pc_type  jacobi lu \
#-ksp_type cg \
#-ksp_atol 1.0e-22 \
#-boundary_2 "X1 001 0 0 1"
#-ksp_type cg      \
#-print_petsc      \
#-part_geom        \
#-normal           \
#-eigen            \
#-ksp_atol 1.0e-10 \
#-ksp_ktol 1.0e-10 \
#-ksp_rtol 1.0e-10 \
