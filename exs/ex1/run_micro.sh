#!/bin/bash


#08-08-2017
function cube_seq {
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
function cube_par {
NM=4

./mpirun -np $NM ../../micro/micro         \
    -input ex1.spu                         \
    -mesh ../../meshes/cube_unif/cube.msh  \
    -mesh_gmsh                             \
    -ksp_type cg                           \
    -ksp_rtol 1.0e-13                      \
    -options_left 0                        \
    -log_trace micro_trace                 \
    -print_vtu
}

#08-08-2017
function barbero_test_seq {
NM=1

./mpirun -np $NM ../../micro/micro            \
    -input ex1.spu                            \
    -mesh ../../meshes/barbero/MESH01/Mesh01  \
    -mesh_alya                                \
    -ksp_type cg                              \
    -ksp_rtol 1.0e-13                         \
    -pc_type  lu                              \
    -options_left 0                           \
    -log_trace micro_trace                    \
    -print_vtu
}

#08-08-2017
function barbero_test_par {
NM=4

./mpirun -np $NM ../../micro/micro            \
    -input ex1.spu                            \
    -mesh ../../meshes/barbero/MESH01/Mesh01  \
    -mesh_alya                                \
    -ksp_type cg                              \
    -ksp_rtol 1.0e-13                         \
    -ksp_atol 1.0e-19                         \
    -options_left 0                           \
    -log_trace micro_trace                    \
    -print_vtu
}

#09-08-2017
function cube_cube_hole_seq {
NM=1

./mpirun -np $NM ../../micro/micro                  \
    -input ex1.spu                                  \
    -mesh ../../meshes/cube_hole/cube_cube_hole.msh \
    -mesh_gmsh                                      \
    -ksp_type cg                                    \
    -ksp_rtol 1.0e-13                               \
    -ksp_atol 1.0e-19                               \
    -pc_type  lu                                    \
    -options_left 0                                 \
    -print_disp                                     \
    -log_trace micro_trace                          \
    -print_vtu
}

#09-08-2017
function cube_cube_hole_par {
NM=4

./mpirun -np $NM ../../micro/micro                  \
    -input ex1.spu                                  \
    -mesh ../../meshes/cube_hole/cube_cube_hole.msh \
    -mesh_gmsh                                      \
    -ksp_type cg                                    \
    -ksp_rtol 1.0e-13                               \
    -ksp_atol 1.0e-19                               \
    -pc_type jacobi                                 \
    -options_left 0                                 \
    -log_trace micro_trace                          \
    -print_vtu
}

#09-08-2017
function cube_chole_fill_seq {
NM=1

./mpirun -np $NM ../../micro/micro                       \
    -input ex1.spu                                       \
    -mesh ../../meshes/cube_hole/cube_cube_hole_fill.msh \
    -mesh_gmsh                                           \
    -pc_type  lu                                         \
    -options_left 0                                      \
    -print_disp                                          \
    -log_trace micro_trace                               \
    -print_vtu
}

#09-08-2017
function cube_chole_fill_par {
NM=4

./mpirun -np $NM ../../micro/micro                       \
    -input ex1.spu                                       \
    -mesh ../../meshes/cube_hole/cube_cube_hole_fill.msh \
    -mesh_gmsh                                           \
    -ksp_type cg                                         \
    -ksp_rtol 1.0e-13                                    \
    -ksp_atol 1.0e-19                                    \
    -pc_type bjacobi                                     \
    -options_left 0                                      \
    -log_trace micro_trace                               \
    -print_vtu
}

#15-08-2017
function cube_basic_par {
NM=4

./mpirun -np $NM ../../micro/micro                       \
    -input ex1.spu                                       \
    -mesh ../../meshes/cube_hole/cube_basic.msh          \
    -mesh_gmsh                                           \
    -ksp_type cg                                         \
    -ksp_rtol 1.0e-13                                    \
    -ksp_atol 1.0e-19                                    \
    -pc_type bjacobi                                     \
    -options_left 0                                      \
    -log_trace micro_trace                               \
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
  -pc_type ilu \
  -pc_factor_nonzeros_along_diagonal \
  -log_trace micro_trace \
  -homo_ld_seq \
  -print_vtu \
  -print_petsc
  #-ksp_type gmres \
  #-ksp_rtol 1.0e-20 \
  #-ksp_atol 1.0e-20 \
  #-pc_type ilu \
}

#cube_seq
#cube_par
#barbero_test_seq
#barbero_test_par
#cube_cube_hole_seq
#cube_cube_hole_par
#cube_chole_fill_seq
#cube_chole_fill_par
#cube_basic_par
#cube_fiber_seq
#cube_fiber_par
#struct_fiber_seq
#struct_fiber_par
#cube_2d_ld_seq
