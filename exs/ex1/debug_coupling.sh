#!/bin/bash

break_mac=( 'mac_main.c:78' ) 
break_mic=( 'mic_main.c:78' ) 
#break_mac=( 'spu_mesh.c:136' ) 
#break_mic=( 'spu_mesh.c:136' ) 
#break_mac=( 'macmic.c:143' ) 
#break_mic=( 'macmic.c:149' ) 
#break_mac=( 'mac_alloc.c:69' ) 
#break_mac=( 'spu_assembly.c:123' ) 
#break_mac=( 'spu_boundary.c:148' ) 

NM=1
Nm=1


# BREAKPOINTS
for i in ${break_mac[@]}
do
  exopt_mac="$exopt_mac -ex 'break $i' "
done
exopt_mac+="-ex 'r'"

for i in ${break_mic[@]}
do
  exopt_mic="$exopt_mic -ex 'break $i' "
done
exopt_mic+="-ex 'r'"

gdbcomm_mac="gdb $exopt_mac --args  ../../macro/macro ex1.spu"
gdbcomm_mac+=" -mesh ../../meshes/cube_unif/cube.msh"
gdbcomm_mac+=" -coupl 1"
gdbcomm_mac+=" -testcomm 1"

gdbcomm_mic="gdb $exopt_mic --args  ../../micro/micro ex1.spu"
gdbcomm_mic+=" -mesh ../../meshes/cube_unif/cube.msh"
gdbcomm_mic+=" -coupl 1"

./mpirun -np $NM xterm -e "$gdbcomm_mac" : -np $Nm xterm -e "$gdbcomm_mic"
