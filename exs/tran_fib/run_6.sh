#!/bin/bash
# 
# 4 simulations of a 2D rve with a fiber transversal to the plane
#
# 1) direct calculation saved on direct/
# 2) multiscale with taylor_s assumtion calculation saved on taylor_s/
# 3) multiscale with taylor_p assumtion calculation saved on taylor_p/
# 4) multiscale with unifst assumtion calculation saved on unifst/
#

if ! [ -d "run_6" ]; then
  mkdir run_6
fi

#---------------------------------------------------------------------

function direct {

NM=4
Nm=1

./mpirun -np $NM ../../macro/macro \
    -input ex2_2.spu \
    -mesh meshes/direct_2/direct.msh \
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

if [ -d "run_6/direct" ]; then
  rm -f run_6/direct/*
else
  mkdir run_6/direct
fi
mv macro_* run_6/direct/.

}

#---------------------------------------------------------------------

function taylor_s {

NM=1
Nm=1

./mpirun \
    -np $NM ../../macro/macro \
    -coupl \
    -input ex2_2.spu \
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
    -input ex2_2.spu \
    -mesh_gmsh \
    -mesh meshes/rve_4/rve_5.msh \
    -dim 2 \
    -pc_type lu \
    -part_geom \
    -homo_taylor_s \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "run_6/taylor_s" ]; then
  rm -f run_6/taylor_s/*
else
  mkdir run_6/taylor_s
fi
mv macro_* run_6/taylor_s/.

}

#---------------------------------------------------------------------

function taylor_p {

NM=1
Nm=1

./mpirun \
    -np $NM ../../macro/macro \
    -coupl \
    -input ex2_2.spu \
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
    -input ex2_2.spu \
    -mesh_gmsh \
    -mesh meshes/rve_4/rve_5.msh \
    -dim 2 \
    -pc_type lu \
    -part_geom \
    -homo_taylor_p \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "run_6/taylor_p" ]; then
  rm -f run_6/taylor_p/*
else
  mkdir run_6/taylor_p
fi
mv macro_* run_6/taylor_p/.

}


#---------------------------------------------------------------------

function unifst {

NM=1
Nm=1

./mpirun \
    -np $NM ../../macro/macro \
    -coupl \
    -input ex2_2.spu \
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
    -input ex2_2.spu \
    -mesh_gmsh \
    -mesh meshes/rve_4/rve_5.msh \
    -dim 2 \
    -pc_type lu \
    -part_geom \
    -homo_us \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "run_6/unifst" ]; then
  rm -f run_6/unifst/*
else
  mkdir run_6/unifst
fi
mv macro_* run_6/unifst/.

}

#---------------------------------------------------------------------

direct
taylor_s
taylor_p
unifst
