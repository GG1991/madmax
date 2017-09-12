#!/bin/bash

NM=1

#break_mac=( 'mac_main.c:xx' ) 
#break_mac=( 'macmic.c:xx' ) 
#break_mac=( 'mac_alloc.c:xx' ) 
#break_mac=( 'spu_mesh.c:xx' ) 
#break_mac=( 'spu_assembly.c:xx' ) 
#break_mac=( 'spu_parser.c:xx' ) 
#break_mac=( 'spu_boundary.c:xx' ) 

# BREAKPOINTS
for i in ${break_mac[@]}
do
  exopt_mac+="$exopt_mac -ex 'break $i' "
done
exopt_mac+="-ex 'r'"

function macro_case {

  NM=1

eval ./mpirun -np $NM xterm -e gdb "$exopt_mac" -q --args ../../macro/macro \
    -input ex1_2d.spu \
    -mesh ../../meshes/cube_fiber/struct_fiber_2d_5_5.msh \
    -dim 2 \
    -mesh_gmsh \
    -ksp_type cg \
    -ksp_rtol 1.0e-13 \
    -log_trace macro_trace \
    -print_vtu \
    -tf 1.0 \
    -dt 1.0 \
    -options_left 0

}

#macro_case
