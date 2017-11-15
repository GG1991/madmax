#!/bin/bash
# 
# 3 simulations
#
# 1) direct calculation saved on direct/
# 2) multiscale with taylor assumtion calculation saved on taylor/
# 3) multiscale with unifst assumtion calculation saved on unifst/
#

if ! [ -d "run_1" ]; then
  mkdir run_1
fi

#---------------------------------------------------------------------

function direct {

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

if [ -d "run_1/direct" ]; then
  rm -f run_1/direct/*
else
  mkdir run_1/direct
fi
mv macro_* run_1/direct/.

}

#---------------------------------------------------------------------

function taylor_s {

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
    -mesh meshes/rve_1/rve_5.msh \
    -dim 2 \
    -pc_type lu \
    -part_geom \
    -homo_taylor_s \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "run_1/taylor_s" ]; then
  rm -f run_1/taylor_s/*
else
  mkdir run_1/taylor_s
fi
mv macro_* run_1/taylor_s/.

}

#---------------------------------------------------------------------

function taylor_p {

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
    -mesh meshes/rve_1/rve_5.msh \
    -dim 2 \
    -pc_type lu \
    -part_geom \
    -homo_taylor_p \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "run_1/taylor_p" ]; then
  rm -f run_1/taylor_p/*
else
  mkdir run_1/taylor_p
fi
mv macro_* run_1/taylor_p/.

}


#---------------------------------------------------------------------

function unifst {

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
    -mesh meshes/rve_1/rve_5.msh \
    -dim 2 \
    -pc_type lu \
    -part_geom \
    -homo_us \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "run_1/unifst" ]; then
  rm -f run_1/unifst/*
else
  mkdir run_1/unifst
fi
mv macro_* run_1/unifst/.

}

#---------------------------------------------------------------------

function unifst_2 {

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
    -mesh meshes/rve_2/rve_1.msh \
    -dim 2 \
    -pc_type lu \
    -part_geom \
    -homo_us \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -fiber_cilin 0.4,0.0,0.0 \
    -options_left 0

if [ -d "run_1/unifst_2" ]; then
  rm -f run_1/unifst_2/*
else
  mkdir run_1/unifst_2
fi
mv macro_* run_1/unifst_2/.

}

#direct
taylor_s
taylor_p
unifst
#unifst_2
