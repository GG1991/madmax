#!/bin/bash
# 
# simulations of a 2D eigenvalue problem
#

if ! [ -d "run_8" ]; then
  mkdir run_8
fi

#---------------------------------------------------------------------

function eigensystem {

NM=1

#xterm -e gdb --args 

./mpirun -np $NM ../../macro/macro \
    -input ex2.spu \
    -dim 2 \
    -mesh_gmsh \
    -mesh meshes/homoge/homog_20.msh \
    -print_vtu \
    -print_petsc \
    -part_geom \
    -eigensys \
    -options_left 0

if [ -d "run_8/direct" ]; then
  rm -f run_8/direct/*
else
  mkdir run_8/direct
fi
mv macro_* run_8/direct/.

}

eigensystem
