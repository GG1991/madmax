#!/bin/bash

NM=1

#break_mic=( 'mic_main.c:354' ) 
#break_mic=( 'spu_mesh.c:xx' ) 
#break_mic=( 'micmic.c:xx' ) 
#break_mic=( 'mic_alloc.c:xx' ) 
#break_mic=( 'spu_assembly.c:52' ) 
#break_mic=( 'mic_boundary.c:xx' ) 

# BREAKPOINTS
for i in ${break_mic[@]}
do
  exopt_mic="$exopt_mac -ex 'break $i' "
done
exopt_mic+="-ex 'r'"

gdbcomm_mic="gdb $exopt_mic --args  ../../micro/micro "
gdbcomm_mic+=" -mesh ../../meshes/cube_unif/cube.msh"

function ex_common {
eval ./mpirun -np $NM xterm -e gdb "$exopt_mic" --args  ../../micro/micro   \
       -input ex1.spu                         \
       -mesh ../../meshes/cube_unif/cube.msh  \
       -ksp_type cg                           \
       -ksp_rtol 1.0e-13                      \
       -options_left 0                        \
       -print_part
}

function ex_valgrind {
eval ./mpirun -np $NM valgrind --log-file=\"valgrind.out\" --leak-check=full ../../micro/micro \
       -input ex1.spu                         \
       -mesh ../../meshes/cube_unif/cube.msh  \
       -ksp_type cg                           \
       -ksp_rtol 1.0e-13                      \
       -options_left 0                        \
       -print_part
}

#ex_common
#ex_valgrind
