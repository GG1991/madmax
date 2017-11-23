#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 

if [ $# -eq 0 ]; then
  NP=1;
else
  NP=$1;
fi

#$MPIEXEC -np $NP xterm -e gdb -x file.gdb --args ../../micro/micro \
$MPIEXEC -np $NP  ../../micro/micro \
    -material      "MATRIX TYPE_0 1.0e7 1.0e6 0.3","FIBER TYPE_0 1.0e7 1.0e7 0.3" \
    -fiber_cilin   0.75,0.0,0.0,0.0   \
    -dim 2                 \
    -struct_n      75,75   \
    -struct_l      3.0,3.0 \
    -pc_type       jacobi  \
    -ksp_type      cg      \
    -print_vtu             \
    -homo_us               \
    -options_left 0

#xterm -e gdb --args 
#xterm -e gdb -x file.gdb --args 
#-pc_type lu \
#-pc_type  jacobi \
#-ksp_type cg \
#-ksp_atol 1.0e-22 \
