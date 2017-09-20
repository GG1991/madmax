#!/bin/bash

NM=1

./mpirun -np $NM ../../macro/macro \
    -input ex2.spu \
    -mesh struct_fiber_2d_10_10.msh \
    -dim 2 \
    -mesh_gmsh \
    -ksp_type cg \
    -ksp_rtol 1.0e-12 \
    -pc_type jacobi \
    -print_vtu \
    -part_geom \
    -nr_norm_tol 1.0e-6 \
    -nr_max_its 2 \
    -tf 1.0 \
    -dt 1.0 \
    -nx_interp 11 \
    -ny_interp 11 \
    -nz_interp 11 \
    -options_left 0

    #-part_meshkway \
    #-ksp_monitor_true_residual \
    #-pc_type lu \
    #-ksp_type cg \
    #-ksp_rtol 1.0e-13 \
    #-log_trace macro_trace \
