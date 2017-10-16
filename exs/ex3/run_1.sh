#!/bin/bash
# 
# simulations of a 2D eigenvalue problem
# a) Simulation of the direct problem
# b) Simulation using multiscale approach
# c) Simulation using parallel mixture theory
# d) Simulation using serial mixture theory
#
# In all the simulation E_i, v_i, v_m remains the same 
# we vary first  E_m from E_m / E_m_0 = 1 to E_m / E_m_0 = 1000
# we vary second r_m from m_m / r_m_0 = 1 to r_m / r_m_0 = 1000
#
# we print the eigenvalue on screen and we take it with awk
#

if ! [ -d "run_1" ]; then
  mkdir run_1
fi

#---------------------------------------------------------------------

function direct {

NM=1

#xterm -e gdb --args 

./mpirun -np $NM ../../macro/macro \
    -input ex2.spu \
    -dim 2 \
    -mesh_gmsh \
    -mesh meshes/direct_1/direct_10.msh \
    -print_vtu \
    -part_geom \
    -eigensys \
    -eps_nev 1 \
    -options_left 0 > macro.out

    #-mesh meshes/homoge/homog_20.msh \
    #-eps_tol 1.0e-7 \
    #-eps_max_it 400 \

if [ -d "run_8/direct" ]; then
  rm -f run_1/direct/*
else
  mkdir run_1/direct
fi
mv macro_* run_1/direct/.

}

eigensystem
