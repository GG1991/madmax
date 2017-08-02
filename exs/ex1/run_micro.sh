#!/bin/bash

NM=4

./mpirun -np $NM ../../micro/micro ex1.spu \
    -mesh ../../meshes/cube_unif/cube.msh  \
    -ksp_type cg                           \
    -ksp_rtol 1.0e-13                      \
    -options_left 0                        \
    -print_part
