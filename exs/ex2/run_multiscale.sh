#!/bin/bash
# 
# Multiscale calculation we test different...
#
#
#

NM=1
Nm=1

function multiscale {

m4 -Dlx_m4=30 -DNx_m4=2 ../../meshes/cube_fiber/struct_homog_2d.geo.m4 > struct_homog_2d.geo
gmsh -2 struct_homog_2d.geo 

./mpirun \
    -np $NM ../../macro/macro \
    -input ex2.spu \
    -mesh_gmsh \
    -mesh struct_homog_2d.msh \
    -dim 2 \
    -pc_type lu \
    -print_vtu \
    -part_geom \
    -nr_norm_tol 1.0e-6 \
    -nr_max_its 2 \
    -tf 1.0 \
    -dt 1.0 \
    -options_left 0
: \
    -np $Nm ../../micro/micro \
    -input ex2.spu \
    -mesh_gmsh \
    -mesh cube_fiber_2d.msh \
    -dim 2 \
    -pc_type lu \
    -print_vtu \
    -part_geom \
    -homo_taylor \
    -nr_norm_tol 1.0e-6 \
    -nr_max_its 2 \
    -options_left 0

    #-part_meshkway \
    #-pc_type lu \
    #-ksp_type cg \
    #-ksp_rtol 1.0e-13 \
    #-log_trace macro_trace \
    #-ksp_monitor_true_residual \

}

multiscale
