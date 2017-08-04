#!/bin/bash

NM=1

./mpirun -np $NM ../../micro/micro         \
    -input ex1.spu                         \
    -mesh ../../meshes/cube_unif/cube.msh  \
    -ksp_type cg                           \
    -ksp_rtol 1.0e-13                      \
    -pc_type  lu                           \
    -options_left 0                        \
    -print_disp                            \
    -print_part

#NM=2
#
#./mpirun -np $NM ../../micro/micro         \
#    -input ex1.spu                         \
#    -mesh ../../meshes/cube_unif/cube.msh  \
#    -ksp_type cg                           \
#    -ksp_rtol 1.0e-13                      \
#    -options_left 0                        \
#    -print_disp                            \
#    -print_part
