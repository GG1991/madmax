#!/bin/bash


break_mic=( 'mic_main.c:144' ) 
#break_mic=( 'spu_mesh.c:1033' ) 
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
       -mesh_gmsh                             \
       -ksp_type cg                           \
       -ksp_rtol 1.0e-13                      \
       -options_left 0                        \
       -print_part
}

function debug_valgrind {
eval ./mpirun -np $NM valgrind --log-file=\"valgrind.out\" --leak-check=full ../../micro/micro \
       -input ex1.spu                         \
       -mesh ../../meshes/cube_unif/cube.msh  \
       -mesh_gmsh                             \
       -ksp_type cg                           \
       -ksp_rtol 1.0e-13                      \
       -options_left 0                        \
       -print_part
}

function barbero_debug {
eval ./mpirun -np $NM xterm -e gdb "$exopt_mic" --args ../../micro/micro   \
    -input ex1.spu                            \
    -mesh ../../meshes/barbero/MESH01/Mesh01  \
    -mesh_gmsh                                \
    -mesh_alya                                \
    -ksp_type cg                              \
    -ksp_rtol 1.0e-13                         \
    -options_left 0                           \
    -print_disp                               \
    -log_trace micro_trace                    \
    -print_part
}

function cube_cube_hole_fill_seq {

NM=1

eval ./mpirun -np $NM xterm -e gdb "$exopt_mic" --args ../../micro/micro   \
    -input ex1.spu                                                         \
    -mesh ../../meshes/cube_hole/cube_cube_hole_fill.msh                   \
    -mesh_gmsh                                                             \
    -pc_type  lu                                                           \
    -options_left 0                                                        \
    -print_disp                                                            \
    -log_trace micro_trace                                                 \
    -print_part
}

#ex_common
#debug_valgrind
#barbero_debug
#cube_cube_hole_fill_seq
