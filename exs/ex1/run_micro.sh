#!/bin/bash


function cube_seq {
NM=1

./mpirun -np $NM ../../micro/micro         \
    -input ex1.spu                         \
    -mesh ../../meshes/cube_unif/cube.msh  \
    -mesh_gmsh                             \
    -ksp_type cg                           \
    -ksp_rtol 1.0e-13                      \
    -pc_type  lu                           \
    -options_left 0                        \
    -print_disp                            \
    -print_part
}


function cube_par {
NM=4

./mpirun -np $NM ../../micro/micro         \
    -input ex1.spu                         \
    -mesh ../../meshes/cube_unif/cube.msh  \
    -mesh_gmsh                             \
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
    -pc_type  lu                              \
    -options_left 0                           \
    -print_disp                               \
    -log_trace macro_trace                    \
    -print_part
}

function barbero_test_par {
NM=4

./mpirun -np $NM ../../micro/micro            \
    -input ex1.spu                            \
    -mesh ../../meshes/barbero/MESH01/Mesh01  \
    -mesh_alya                                \
    -ksp_type cg                              \
    -ksp_rtol 1.0e-13                         \
    -ksp_atol 1.0e-19                         \
    -options_left 0                           \
    -print_disp                               \
    -log_trace macro_trace                    \
    -print_part
}

#cube_seq
#cube_par
#barbero_test_seq
#barbero_test_par
