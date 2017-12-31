#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 

#$MPIEXEC -np 1 xterm -e gdb -x file.gdb --args ../../micro/micro \
$MPIEXEC -np 1  ../../micro/micro \
    -struct_n 75,75 \
    -dim 2 \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3" \
    -micro_struct "fiber_cilin 3.0 3.0 1 1 0.75 0.0 0.0" \
    -pc_type jacobi \
    -ksp_type cg \
    -homo_us \
    -print_pvtu

#xterm -e gdb --args 
#xterm -e gdb -x file.gdb --args 
#-pc_type lu \
#-pc_type  jacobi \
#-ksp_type cg \
#-ksp_atol 1.0e-22 \
#
#-micro_struct  "size"[dim]     \
#               "nx_fib"        \
#               "ny_fib"        \
#               "radio"         \
#               "desv"[2]
#
#-micro_struct  "fiber_line     \
#               "size"[dim]     \
#               "ntype"         \
#               "nfib[ntype]"   \
#               "tetha[ntype]"  \
#               "seps[ntype]"   \
#               "width[ntype]"  \
#               "desv[ntype]"

#-micro_struct  "fiber_cilin 3.0 3.0 1 1 0.75 0.0 0.0" \
#-micro_struct  "fiber_line 3.0 3.0 1 1 0.0 0.4 1.0 0.0" \
#-micro_struct  "fiber_line 3.0 3.0 2 9 9 0.785398 2.35619 0.4 0.4 0.2 0.2 0.0 0.0" \
