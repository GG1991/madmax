#!/bin/bash


function struct_fiber_2d {

NM=4
./mpirun -np $NM ../../macro/macro \
    -input ex2.spu \
    -mesh struct_fiber_2d_10_10.msh \
    -dim 2 \
    -mesh_gmsh \
    -ksp_type bicg \
    -ksp_rtol 1.0e-8 \
    -ksp_monitor_true_residual \
    -pc_type jacobi \
    -print_vtu \
    -part_geom \
    -nr_norm_tol 1.0e-6 \
    -nr_max_its 2 \
    -tf 2.0 \
    -dt 1.0 \
    -options_left 0

    #-part_meshkway \
    #-pc_type lu \
    #-ksp_type cg \
    #-ksp_rtol 1.0e-13 \
    #-log_trace macro_trace \

}

function debug_struct_fiber_2d {

NM=1
./mpirun -np $NM xterm -e gdb "$exopt_mic" --args ../../macro/macro \
    -input ex2.spu \
    -mesh struct_fiber_2d_10_10.msh \
    -dim 2 \
    -mesh_gmsh \
    -ksp_type cg \
    -ksp_rtol 1.0e-13 \
    -print_vtu \
    -nr_norm_tol 1.0e-4 \
    -tf 2.0 \
    -dt 1.0 \
    -options_left 0

    #-pc_type lu \
    #-ksp_type cg \
    #-ksp_rtol 1.0e-13 \
    #-log_trace macro_trace \

}

#struct_fiber_2d
#debug_struct_fiber_2d
