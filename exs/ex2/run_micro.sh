#!/bin/bash

NM=1

./mpirun -np $NM ../../micro/micro \
    -input ex2.spu \
    -mesh cube_fiber_2d.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_unif_strains \
    -print_vtu

#./mpirun -np $NM xterm -e gdb --args ../../micro/micro \
