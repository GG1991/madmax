#!/bin/bash
# 
# Multiscale calculation we test different...
#
#
#

NM=1
Nm=1

m4 -Dlx_m4=30 -DNx_m4=Nx_m4_m4 ../../meshes/cube_fiber/struct_homog_2d.geo.m4 > struct_homog_2d.geo
gmsh -2 struct_homog_2d.geo > /tmp/null

./mpirun \
    -np $NM ../../macro/macro \
    -coupl \
    -input ex2.spu \
    -mesh_gmsh \
    -mesh struct_homog_2d.msh \
    -dim 2 \
    -pc_type lu \
    -print_vtu \
    -part_geom \
    -nr_norm_tol 1.0e-6 \
    -nr_max_its 3 \
    -tf 1.0 \
    -dt 1.0 \
    -options_left 0 \
: \
    -np $Nm ../../micro/micro \
    -coupl \
    -input ex2.spu \
    -mesh_gmsh \
    -mesh cube_fiber_2d.msh \
    -dim 2 \
    -pc_type lu \
    -part_geom \
    -homo_m4 \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0
#    #xterm -e gdb --args
#    #-part_meshkway \
#    #-pc_type lu \
#    #-ksp_type cg \
#    #-ksp_rtol 1.0e-13 \
#    #-log_trace macro_trace \
#    #-ksp_monitor_true_residual \
