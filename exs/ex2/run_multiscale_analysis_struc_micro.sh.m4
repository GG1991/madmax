#!/bin/bash
# 
# Multiscale calculation we test different...
#
#
#

NM=1
Nm=1

m4 -Dlx_m4=30 -DNx_m4=21 ../../meshes/cube_fiber/struct_homog_2d.geo.m4 > struct_homog_2d.geo
gmsh -2 struct_homog_2d.geo > /tmp/null
m4 -Dlc_m4=lc_m4_m4 -DN_m4=N_m4_m4 ../../meshes/cube_unif/cube_2d_analysis.geo.m4 > cube_2d.geo
gmsh -2 cube_2d.geo > /tmp/null

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
    -mesh cube_2d.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_unif_strains \
    -print_vtu \
    -fiber_cilin 0.4,-0.0,+0.0 \
    -fiber_nx nx_m4 \
    -fiber_ny nx_m4 \
    -options_left 0

#    #xterm -e gdb --args
#    #-part_meshkway \
#    #-pc_type lu \
#    #-ksp_type cg \
#    #-ksp_rtol 1.0e-13 \
#    #-log_trace macro_trace \
#    #-ksp_monitor_true_residual \
