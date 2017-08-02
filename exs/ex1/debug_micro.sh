#!/bin/bash

NM=1

#break_mic=( 'mic_main.c:281' ) 
#break_mic=( 'spu_mesh.c:136' ) 
#break_mic=( 'micmic.c:143' ) 
#break_mic=( 'mic_alloc.c:69' ) 
#break_mic=( 'spu_assembly.c:123' ) 
break_mic=( 'mic_boundary.c:373' ) 

# BREAKPOINTS
for i in ${break_mic[@]}
do
  exopt_mic="$exopt_mac -ex 'break $i' "
done
exopt_mic+="-ex 'r'"

gdbcomm_mic="gdb $exopt_mic --args  ../../micro/micro "
gdbcomm_mic+=" -mesh ../../meshes/cube_unif/cube.msh"
eval ./mpirun -np $NM xterm -e gdb "$exopt_mic" --args  ../../micro/micro   \
       -input ex1.spu                         \
       -mesh ../../meshes/cube_unif/cube.msh  \
       -ksp_type cg                           \
       -ksp_rtol 1.0e-13                      \
       -options_left 0                        \
       -print_part

#eval ./mpirun -np $NM valgrind --log-file=\"valgrind.out\" --leak-check=full ../../micro/micro \
#       -input ex1.spu                         \
#       -mesh ../../meshes/cube_unif/cube.msh  \
#       -ksp_type cg                           \
#       -ksp_rtol 1.0e-13                      \
#       -options_left 0                        \
#       -print_part
#
