#!/bin/bash
# 
# simulations of a 2D eigenvalue problem
# a) Simulation of the direct problem
# b) Simulation using multiscale approach
# c) Simulation using parallel mixture theory
# d) Simulation using serial mixture theory
#
# In all the simulation E_i, v_i, rho_i, v_m, rho_m remains the same 
# we vary E_m from E_m / E_m_0 = 1 to E_m / E_m_0 = 1000
#
# we print the eigenvalue on screen and we take it with awk
#

if ! [ -d "run_1" ]; then
  mkdir run_1
fi

E_m=( 2.0e5 5.0e5 8.0e5 1.0e6 2.0e6 5.0e6 8.0e6 1.0e7 2.0e7 5.0e7 8.0e7 1.0e8 \
2.0e8 5.0e8 )
nev=3 # number of request eigenvalues

#---------------------------------------------------------------------

function direct {

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
      -eps_nev $nev \
      -options_left 0 > macro.out
      #-mesh meshes/homoge/homog_20.msh \
      #-eps_tol 1.0e-7 \
      #-eps_max_it 400 \
      #-eps_view \

  if [ -d "run_1/direct_$i" ]; then
    rm -f run_1/direct_$i/*
  else
    mkdir run_1/direct_$i
  fi
  mv macro* ex.spu run_1/direct_$i/.
  echo "run_1/direct_$i done"

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
      -eps_nev $nev \
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
      -eps_nev $nev \
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
      -eps_nev $nev \
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

#---------------------------------------------------------------------

function extract_data {

# extract results 

rm -f em.dat
for i in `seq 1 ${#E_m[@]}`; do
 echo ${E_m[$((i-1))]} >> em.dat
done
cp em.dat omega_aux.dat

for v in `seq 1 $nev`; do

 rm -f omega_tp_$v.dat omega_ts_$v.dat omega_direct_$v.dat omega_us_$v.dat 

 for i in `seq 1 ${#E_m[@]}`; do
  sed -n "/omega $((v-1))/p" run_1/direct_$i/macro.out   | awk '{print $4}' >> omega_direct_$v.dat
 done
 for i in `seq 1 ${#E_m[@]}`; do
  sed -n "/omega $((v-1))/p" run_1/multi_us_$i/macro.out | awk '{print $4}' >> omega_us_$v.dat
 done
 for i in `seq 1 ${#E_m[@]}`; do
  sed -n "/omega $((v-1))/p" run_1/multi_tp_$i/macro.out | awk '{print $4}' >> omega_tp_$v.dat
 done
 for i in `seq 1 ${#E_m[@]}`; do
  sed -n "/omega $((v-1))/p" run_1/multi_ts_$i/macro.out | awk '{print $4}' >> omega_ts_$v.dat
 done

 paste omega_aux.dat omega_direct_$v.dat omega_us_$v.dat omega_tp_$v.dat omega_ts_$v.dat \
 > omega_aux_1.dat ; mv omega_aux_1.dat omega_aux.dat
 
 mv omega_direct_$v.dat omega_us_$v.dat omega_tp_$v.dat omega_ts_$v.dat run_1/.

done
mv omega_aux.dat run_1/omega_vs_em.dat
mv em.dat run_1/.

}

#---------------------------------------------------------------------

direct
multi_us
multi_tp
multi_ts
extract_data
