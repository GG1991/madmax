#!/bin/bash

NM=1

# 16-08-2017
# Cube with cilindrical fiber in the middle
function cube_fiber_2d {

./mpirun -np $NM ../../micro/micro \
    -input ex2.spu \
    -mesh cube_fiber_2d.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_ld \
    -print_vtu
}

cube_fiber_2d
