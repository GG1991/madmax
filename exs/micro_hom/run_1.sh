#!/bin/bash

#./mpiexec -np 1 xterm -e gdb -x file.gdb --args  ../../micro/micro \
./mpiexec -np 1 ../../micro/micro \
    -dim 2 \
    -struct_n    15,15 \
    -struct_l    0.1,0.1 \
    -fiber_cilin 0.4,0.0,0.0,0.0   \
    -mat_fiber_t0  1.0e6,1.0e7,0.3 \
    -mat_matrix_t0 1.0e6,1.0e6,0.3 \
    -pc_type lu \
    -print_vtu \
    -homo_us \
    -options_left 0
#: -np 1 xterm -e gdb -x file.gdb --args ../../micro/micro \
#   -dim 2 \
#   -struct_n    5,5 \
#   -struct_l    0.1,0.1 \
#   -fiber_cilin 0.4,0.0,0.0,0.0   \
#   -mat_fiber_t0  1.0e6,1.0e7,0.3 \
#   -mat_matrix_t0 1.0e6,1.0e6,0.3 \
#   -pc_type lu \
#   -print_vtu \
#   -homo_us \
#   -options_left 0

#xterm -e gdb --args 
#xterm -e gdb -x file.gdb --args 
