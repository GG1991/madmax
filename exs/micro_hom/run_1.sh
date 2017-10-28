#!/bin/bash

./mpiexec -np 1 ../../micro/micro \
    -input ex2.spu \
    -mesh_gmsh \
    -mesh meshes/rve_struc/rve_1.msh \
    -fiber_cilin 0.4,0.0,0.0,0.0   \
    -mat_fiber_t0  1.0e6,1.0e7,0.3 \
    -mat_matrix_t0 1.0e6,1.0e6,0.3 \
    -dim 2 \
    -pc_type lu \
    -print_vtu \
    -homo_us \
    -options_left 0

#xterm -e gdb --args 
