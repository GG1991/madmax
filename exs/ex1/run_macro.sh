#!/bin/bash

NM=2

./mpirun -np $NM ../../macro/macro ex1.spu \
    -mesh ../../meshes/cube_unif/cube.msh  \
    -ksp_type cg                           \
    -ksp_rtol 1.0e-13                      \
    -options_left 0                        \
    -print_part
