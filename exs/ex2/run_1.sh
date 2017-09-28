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

function taylor_1 {

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
    -homo_taylor_1 \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "run_1/taylor_1" ]; then
  rm -f run_1/taylor_1/*
else
  mkdir run_1/taylor_1
fi
mv macro_* run_1/taylor_1/.

}

#---------------------------------------------------------------------

function taylor_2 {

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
    -homo_taylor_2 \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "run_1/taylor_2" ]; then
  rm -f run_1/taylor_2/*
else
  mkdir run_1/taylor_2
fi
mv macro_* run_1/taylor_2/.

}


#---------------------------------------------------------------------

function unifst_1 {

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
    -homo_unif_strains \
    -nr_norm_tol 1.0e-8 \
    -nr_max_its 3 \
    -options_left 0

if [ -d "run_1/unifst_1" ]; then
  rm -f run_1/unifst_1/*
else
  mkdir run_1/unifst_1
fi
mv macro_* run_1/unifst_1/.

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
    -homo_unif_strains \
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
taylor_1
taylor_2
unifst_1
unifst_2
