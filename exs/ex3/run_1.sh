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

#---------------------------------------------------------------------

function multi_us {

NM=1
#xterm -e gdb --args 
for i in `seq 1 ${#E_m[@]}`; do

  m4 -Drho_m=1.0e6 -DE_m=${E_m[$((i-1))]} ex.spu.m4 > ex.spu
  ./mpirun -np $NM ../../macro/macro \
      -coupl \
      -input ex.spu \
      -dim 2 \
      -mesh_gmsh \
      -mesh meshes/homoge/homog_10.msh \
      -print_vtu \
      -part_geom \
      -eigensys \
      -eps_nev 1 \
      -options_left 0 \
   : -np $NM ../../micro/micro \
      -coupl \
      -input ex.spu \
      -dim 2 \
      -mesh_gmsh \
      -mesh meshes/rve_1/rve_5.msh \
      -part_geom \
      -homo_us \
      -options_left 0 > macro.out
      #-mesh meshes/homoge/homog_20.msh \
      #-eps_tol 1.0e-7 \
      #-eps_max_it 400 \

  if [ -d "run_1/multi_us_$i" ]; then
    rm -f run_1/multi_us_$i/*
  else
    mkdir run_1/multi_us_$i
  fi
  mv macro* ex.spu run_1/multi_us_$i/.
  echo "run_1/multi_us_$i done"

done
}

#---------------------------------------------------------------------

function multi_tp {

NM=1
#xterm -e gdb --args 
for i in `seq 1 ${#E_m[@]}`; do

  m4 -Drho_m=1.0e6 -DE_m=${E_m[$((i-1))]} ex.spu.m4 > ex.spu
  ./mpirun -np $NM ../../macro/macro \
      -coupl \
      -input ex.spu \
      -dim 2 \
      -mesh_gmsh \
      -mesh meshes/homoge/homog_10.msh \
      -print_vtu \
      -part_geom \
      -eigensys \
      -eps_nev 1 \
      -options_left 0 \
   : -np $NM ../../micro/micro \
      -coupl \
      -input ex.spu \
      -dim 2 \
      -mesh_gmsh \
      -mesh meshes/rve_1/rve_5.msh \
      -part_geom \
      -homo_taylor_p \
      -options_left 0 > macro.out
      #-mesh meshes/homoge/homog_20.msh \
      #-eps_tol 1.0e-7 \
      #-eps_max_it 400 \

  if [ -d "run_1/multi_tp_$i" ]; then
    rm -f run_1/multi_tp_$i/*
  else
    mkdir run_1/multi_tp_$i
  fi
  mv macro* ex.spu run_1/multi_tp_$i/.
  echo "run_1/multi_tp_$i done"

done
}

#---------------------------------------------------------------------

function multi_ts {

NM=1
#xterm -e gdb --args 
for i in `seq 1 ${#E_m[@]}`; do

  m4 -Drho_m=1.0e6 -DE_m=${E_m[$((i-1))]} ex.spu.m4 > ex.spu
  ./mpirun -np $NM ../../macro/macro \
      -coupl \
      -input ex.spu \
      -dim 2 \
      -mesh_gmsh \
      -mesh meshes/homoge/homog_10.msh \
      -print_vtu \
      -part_geom \
      -eigensys \
      -eps_nev 1 \
      -options_left 0 \
   : -np $NM ../../micro/micro \
      -coupl \
      -input ex.spu \
      -dim 2 \
      -mesh_gmsh \
      -mesh meshes/rve_1/rve_5.msh \
      -part_geom \
      -homo_taylor_s \
      -options_left 0 > macro.out
      #-mesh meshes/homoge/homog_20.msh \
      #-eps_tol 1.0e-7 \
      #-eps_max_it 400 \

  if [ -d "run_1/multi_ts_$i" ]; then
    rm -f run_1/multi_ts_$i/*
  else
    mkdir run_1/multi_ts_$i
  fi
  mv macro* ex.spu run_1/multi_ts_$i/.
  echo "run_1/multi_ts_$i done"

done
}

function extract_data {

# extract results 
rm -f omega_direct.dat omega_us.dat em.dat run_1/omega_vs_em.dat \
omega_tp.dat omega_ts.dat

for i in `seq 1 ${#E_m[@]}`; do
 echo ${E_m[$((i-1))]} >> em.dat
done

for i in `seq 1 ${#E_m[@]}`; do
 awk '/omega/{print $4}' run_1/direct_1_$i/macro.out >> omega_direct.dat
done

for i in `seq 1 ${#E_m[@]}`; do
 awk '/omega/{print $4}' run_1/multi_us_$i/macro.out >> omega_us.dat
done

for i in `seq 1 ${#E_m[@]}`; do
 awk '/omega/{print $4}' run_1/multi_tp_$i/macro.out >> omega_tp.dat
done

for i in `seq 1 ${#E_m[@]}`; do
 awk '/omega/{print $4}' run_1/multi_ts_$i/macro.out >> omega_ts.dat
done

paste em.dat          omega_direct.dat   > omega_aux_1.dat
mv omega_aux_1.dat    omega_aux.dat
paste omega_aux.dat   omega_us.dat       > omega_aux_1.dat
mv omega_aux_1.dat    omega_aux.dat
paste omega_aux.dat   omega_tp.dat       > omega_aux_1.dat
mv omega_aux_1.dat    omega_aux.dat
paste omega_aux.dat   omega_ts.dat       > omega_aux_1.dat

mv omega_aux_1.dat run_1/omega_vs_em.dat

}

#direct_1
#multi_us
#multi_tp
#multi_ts
extract_data
