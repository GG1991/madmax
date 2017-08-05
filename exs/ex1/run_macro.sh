#!/bin/bash


function nproc_4 {

NM=4
./mpirun -np $NM ../../macro/macro         \
    -input ex1.spu                         \
    -mesh ../../meshes/cube_unif/cube.msh  \
    -ksp_type cg                           \
    -ksp_rtol 1.0e-13                      \
    -options_left 0                        \
    -log_trace macro_trace                 \
    -print_part
}
nproc_4
