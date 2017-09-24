#!/bin/bash
# 
# 3 simulations
#
# 1) direct calculation saved on direct/
# 2) multiscale with taylor assumtion calculation saved on taylor/
# 3) multiscale with unifst assumtion calculation saved on unifst/
#


#---------------------------------------------------------------------

NM=4
Nm=1

./mpirun -np $NM ../../macro/macro \
    -input ex2.spu \
    -mesh meshes/direct/direct.msh \
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

if [ -d "direct" ]; then
  rm -f direct/*
else
  mkdir direct
fi
mv macro_* direct/.

#---------------------------------------------------------------------

NM=1
Nm=1

./mpirun \
    -np $NM ../../macro/macro \
    -coupl \
    -input ex2.spu \
    -mesh_gmsh \
    -mesh meshes/homoge/homog_12.msh \
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
    -mesh meshes/rve/rve_5.msh \
    -dim 2 \
    -pc_type lu \
    -part_geom \
    -homo_taylor \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "taylor" ]; then
  rm -f taylor/*
else
  mkdir taylor
fi
mv macro_* taylor/.

#---------------------------------------------------------------------

NM=1
Nm=1

./mpirun \
    -np $NM ../../macro/macro \
    -coupl \
    -input ex2.spu \
    -mesh_gmsh \
    -mesh meshes/homoge/homog_12.msh \
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
    -mesh meshes/rve/rve_5.msh \
    -dim 2 \
    -pc_type lu \
    -part_geom \
    -homo_unif_strains \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "unifst" ]; then
  rm -f unifst/*
else
  mkdir unifst
fi
mv macro_* unifst/.


