#!/bin/bash


function cube {

NM=4
./mpirun -np $NM ../../macro/macro         \
    -input ex1.spu                         \
    -mesh ../../meshes/cube_unif/cube.msh  \
    -ksp_type cg \
    -ksp_rtol 1.0e-13 \
    -log_trace macro_trace \
    -print_pvtui \
    -options_left 0
}

function struct_fiber {

NM=4
./mpirun -np $NM ../../macro/macro \
    -input ex1.spu \
    -mesh ../../meshes/cube_fiber/struct_fiber.msh \
    -mesh_gmsh \
    -ksp_type cg \
    -ksp_rtol 1.0e-13 \
    -log_trace macro_trace                 \
    -print_vtu \
    -tf 1.0 \
    -dt 1.0 \
    -options_left 0
}
#struct_fiber
