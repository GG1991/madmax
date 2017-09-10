#!/bin/bash

#08-08-2017
function cube {

NM=1

./mpirun -np $NM ../../micro/micro         \
    -input ex1.spu                         \
    -mesh ../../meshes/cube_unif/cube.msh  \
    -mesh_gmsh                             \
    -ksp_type cg                           \
    -ksp_rtol 1.0e-13                      \
    -pc_type  lu                           \
    -options_left 0                        \
    -print_vtu
}

#08-08-2017
function barbero {
NM=1

./mpirun -np $NM ../../micro/micro            \
    -input ex1.spu                            \
    -mesh ../../meshes/barbero/MESH01/Mesh01  \
    -mesh_alya                                \
    -ksp_type cg                              \
    -ksp_rtol 1.0e-13                         \
    -pc_type  lu                              \
    -options_left 0                           \
    -homo_ld \
    -log_trace micro_trace                    \
    -print_vtu
}


#15-08-2017
function cube_hole {

NM=1

./mpirun -np $NM ../../micro/micro \
    -input ex1_2d.spu \
    -mesh ../../meshes/cube_hole/cube_hole.msh \
    -dim 3 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_ld \
    -log_trace micro_trace \
    -print_vtu
}

#15-08-2017
function cube_hole_2d {

NM=1

./mpirun -np $NM ../../micro/micro \
    -input ex1_2d.spu \
    -mesh ../../meshes/cube_hole/cube_hole_2d.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_ld \
    -log_trace micro_trace \
    -print_vtu
}

# 16-08-2017
# Cube with cilindrical fiber in the middle
function cube_fiber_seq {
NM=1

./mpirun -np $NM ../../micro/micro \
    -input ex1.spu \
    -mesh ../../meshes/cube_fiber/cube_fiber.msh \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -log_trace micro_trace \
    -homo_exp \
    -print_vtu
}

# 16-08-2017
# Cube with cilindrical fiber in the middle
function cube_fiber_par {

  NM=4

    ./mpirun -np $NM ../../micro/micro \
    -input ex1.spu \
    -mesh ../../meshes/cube_fiber/cube_fiber.msh \
    -mesh_gmsh \
    -ksp_type cg \
    -ksp_rtol 1.0e-13 \
    -ksp_atol 1.0e-19 \
    -pc_type bjacobi \
    -options_left 0 \
    -log_trace micro_trace \
    -homo_linear \
    -print_vtu
}

function struct_fiber_par {

  NM=4

    ./mpirun -np $NM ../../micro/micro \
    -input ex1.spu \
    -mesh ../../meshes/cube_fiber/struct_fiber.msh \
    -mesh_gmsh \
    -ksp_type cg \
    -ksp_rtol 1.0e-15 \
    -ksp_atol 1.0e-11 \
    -pc_type bjacobi  \
    -options_left 0 \
    -log_trace micro_trace \
    -homo_linear \
    -print_vtu
}

function cube_fiber_ld_seq {
NM=1 # this technique is sequencial only 

  ./mpirun -np $NM ../../micro/micro \
  -input ex1.spu \
  -mesh ../../meshes/cube_fiber/cube_fiber.msh \
  -mesh_gmsh \
  -ksp_type cg \
  -ksp_rtol 1.0e-13 \
  -ksp_atol 1.0e-19 \
  -pc_type bjacobi  \
  -options_left 0 \
  -log_trace micro_trace \
  -homo_ld_seq \
  -print_vtu
}

# 07-09-2017
# Cube 2D
# Linear Displacements with Lagrangian BC
# this technique is sequencial only
function cube_2d_ld_seq {
NM=1  

  ./mpirun -np $NM ../../micro/micro \
  -input ex1_2d.spu \
  -mesh ../../meshes/cube_unif/cube_2d.msh \
  -dim 2 \
  -mesh_gmsh \
  -options_left 0 \
  -pc_type lu \
  -log_trace micro_trace \
  -homo_ld \
  -print_vtu \
  -print_petsc
  #-ksp_type gmres \
  #-ksp_rtol 1.0e-20 \
  #-ksp_atol 1.0e-20 \
  #-pc_type ilu \
  #-homo_ld_seq \
  #-pc_factor_nonzeros_along_diagonal \
}

# 08-09-2017
# Cube 2D with fiber
# Linear Displacements with Lagrangian BC
# this technique is sequencial only
function cube_fiber_2d_ld_seq {
NM=1  

  ./mpirun -np $NM ../../micro/micro \
  -input ex1_2d.spu \
  -mesh ../../meshes/cube_fiber/cube_fiber_2d.msh \
  -dim 2 \
  -mesh_gmsh \
  -options_left 0 \
  -pc_type lu \
  -log_trace micro_trace \
  -homo_ld \
  -print_vtu \
  -print_petsc
  #-ksp_type gmres \
  #-ksp_rtol 1.0e-20 \
  #-ksp_atol 1.0e-20 \
  #-pc_type ilu \
  #-homo_ld_l_seq \
  #-pc_factor_nonzeros_along_diagonal \
}

#cube
#barbero
cube_hole_2d
#cube_hole
