#!/bin/bash


function nproc_1 {
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
}


function nproc_2 {
NM=2

./mpirun -np $NM ../../micro/micro         \
    -input ex1.spu                         \
    -mesh ../../meshes/cube_unif/cube.msh  \
    -ksp_type cg                           \
    -ksp_rtol 1.0e-13                      \
    -options_left 0                        \
    -print_disp                            \
    -log_trace macro_trace                 \
    -print_part
}

function barbero_test_seq {
NM=1

./mpirun -np $NM ../../micro/micro            \
    -input ex1.spu                            \
    -mesh ../../meshes/barbero/MESH01/Mesh01  \
    -mesh_alya                                \
    -ksp_type cg                              \
    -ksp_rtol 1.0e-13                         \
    -options_left 0                           \
    -print_disp                               \
    -log_trace macro_trace                    \
    -print_part
}

#nproc_1
#nproc_2
#barbero_test_seq
