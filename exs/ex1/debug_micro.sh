#!/bin/bash

#break_mic=( 'mic_main.c:xx' ) 
#break_mic=( 'spu_mesh.c:xx' ) 
#break_mic=( 'micmic.c:xx' ) 
#break_mic=( 'mic_alloc.c:xx' ) 
#break_mic=( 'spu_assembly.c:xx' ) 
#break_mic=( 'mic_boundary.c:xx' ) 

# BREAKPOINTS
exopt_mic=''
for i in ${break_mic[@]}
do
  exopt_mic+="$exopt_mac -ex 'break $i' "
done
exopt_mic+="-ex 'r'"

function valgrind {
eval ./mpirun -np $NM valgrind --log-file=\"valgrind.out\" --leak-check=full ../../micro/micro \
       -input ex1.spu                         \
       -mesh ../../meshes/cube_unif/cube.msh  \
       -mesh_gmsh                             \
       -ksp_type cg                           \
       -ksp_rtol 1.0e-13                      \
       -options_left 0                        \
       -print_part
}

function barbero {
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

# 23-08-2017
# Cube with cilindrical hole in the middle
function cube_2d {
NM=1

eval ./mpirun -np $NM xterm -e gdb "$exopt_mic" -q --args ../../micro/micro \
  -input ex1_2d.spu \
  -mesh ../../meshes/cube_unif/cube_2d.msh \
  -dim 2 \
  -mesh_gmsh \
  -pc_type lu \
  -options_left 0 \
  -log_trace micro_trace \
  -homo_taylor \
  -fiber_cilin 0.4,-0.0,+0.0 \
  -fiber_nx 2 \
  -fiber_ny 2 \
  -print_vtu
}

# 23-08-2017
# Cube with cilindrical fiber in the middle
function cube_fiber {
NM=1

eval ./mpirun -np $NM xterm -e gdb "$exopt_mic" -q --args ../../micro/micro \
  -input ex1_2d.spu \
  -mesh ../../meshes/cube_fiber/cube_fiber_2d.msh \
  -dim 2 \
  -mesh_gmsh \
  -pc_type lu \
  -options_left 0 \
  -log_trace micro_trace \
  -homo_ld \
  -print_vtu
}

#valgrind
#barbero
#cube_2d
#cube_fiber
