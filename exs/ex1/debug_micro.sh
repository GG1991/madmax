#!/bin/bash

NM=1

#break_mic=( 'mic_main.c:200' ) 
#break_mic=( 'spu_mesh.c:136' ) 
#break_mic=( 'micmic.c:143' ) 
#break_mic=( 'mic_alloc.c:69' ) 
#break_mic=( 'spu_assembly.c:123' ) 
break_mic=( 'mic_boundary.c:259' ) 

# BREAKPOINTS
for i in ${break_mic[@]}
do
  exopt_mic="$exopt_mac -ex 'break $i' "
done
exopt_mic+="-ex 'r'"

gdbcomm_mic="gdb $exopt_mic --args  ../../micro/micro ex1.spu "
gdbcomm_mic+=" -mesh ../../meshes/cube_unif/cube.msh"
./mpirun -np $NM xterm -e "$gdbcomm_mic"
