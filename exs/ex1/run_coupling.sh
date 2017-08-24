#!/bin/bash
#
#      -ksp_atol <abstol> 	- Sets abstol
#      -ksp_rtol <rtol> 	- Sets rtol
#      -ksp_divtol <dtol> 	- Sets dtol
#      -ksp_max_it <maxits> 	- Sets maxits 

NM=1
Nm=1

function cube_fiber {

./mpirun -np $NM ../../macro/macro \
    -input ex1.spu \
    -mesh ../../meshes/cube_fiber/struct_fiber_homog.msh \
    -mesh_gmsh \
    -ksp_type cg \
    -ksp_rtol 1.0e-13 \
    -ksp_atol 1.0e-11 \
    -coupl \
    -tf 1.0 \
    -dt 1.0 \
    -print_vtu \
    -options_left 0 \
: -np $Nm ../../micro/micro \
    -input ex1.spu \
    -mesh ../../meshes/cube_fiber/cube_fiber.msh \
    -mesh_gmsh \
    -coupl \
    -homo_linear \
    -options_left 0
}

cube_fiber
