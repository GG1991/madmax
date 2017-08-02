#!/bin/bash

NM=1

break_mac=( 'mac_main.c:54' ) 
#break_mac=( 'spu_mesh.c:136' ) 
#break_mac=( 'macmic.c:143' ) 
#break_mac=( 'mac_alloc.c:69' ) 
#break_mac=( 'spu_assembly.c:123' ) 
#break_mac=( 'spu_boundary.c:148' ) 

# BREAKPOINTS
for i in ${break_mac[@]}
do
  exopt_mac="$exopt_mac -ex 'break $i' "
done
exopt_mac+="-ex 'r'"

gdbcomm_mac="gdb $exopt_mac --args  ../../macro/macro ex1.spu -coupl"
gdbcomm_mac+=" -mesh ../../meshes/cube_unif/cube.msh"
./mpirun -np $NM xterm -e "$gdbcomm_mac"
