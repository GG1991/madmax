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

E_m=( 1.0e7 2.0e7 5.0e7 8.0e7 1.0e8 2.0e8 5.0e8 8.0e8 1.0e9 )

#---------------------------------------------------------------------

function direct_1 {

NM=1

#xterm -e gdb --args 


for i in `seq 1 ${#E_m[@]}`; do

  m4 -Drho_m=1.0e6 -DE_m=${E_m[$((i-1))]} ex.spu.m4 > ex.spu
  ./mpirun -np $NM ../../macro/macro \
      -input ex.spu \
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

  if [ -d "run_1/direct_1_$i" ]; then
    rm -f run_1/direct_1_$i/*
  else
    mkdir run_1/direct_1_$i
  fi
  mv macro* ex.spu run_1/direct_1_$i/.
  echo "run_1/direct_1_$i done"

done
}

function ext_direct_1 {

# extract results 
rm -f omega.dat em.dat run_1/omega_vs_em.dat

for i in `seq 1 ${#E_m[@]}`; do
 awk '/omega/{print $4}' run_1/direct_1_$i/macro.out >> omega.dat
 echo ${E_m[$((i-1))]} >> em.dat
done
paste em.dat omega.dat > run_1/omega_vs_em.dat

}
direct_1
ext_direct_1
