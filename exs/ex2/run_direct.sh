#!/bin/bash


function struct_fiber_2d {

NM=4
./mpirun -np $NM ../../macro/macro \
    -input ex2.spu \
    -mesh ../../meshes/cube_fiber/struct_fiber_2d_10_10.msh \
    -dim 2 \
    -mesh_gmsh \
    -ksp_type cg \
    -ksp_rtol 1.0e-13 \
    -log_trace macro_trace \
    -print_vtu \
    -tf 2.0 \
    -dt 1.0 \
    -options_left 0
}

struct_fiber_2d
