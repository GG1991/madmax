#!/bin/bash
#
# Performs an homogenization for different mesh and prints
# the homogenized constitutive tensor
#

NM=1

if [ -d "run_4" ]; then
   rm -rf run_4/*
else
   mkdir run_4
fi

if [ -e "run_4_a.dat" ]; then
   rm run_4_a.dat
fi
if [ -e "run_4_b.dat" ]; then
   rm run_4_b.dat
fi
if [ -e "run_4_c.dat" ]; then
   rm run_4_c.dat
fi

#---------------------------------------------------------------------- 

for i in $(seq 1 5);
do

    echo "case $i"
    ./mpirun -np $NM ../../micro/micro \
    -input ex2.spu \
    -mesh meshes/rve_3/rve_$i.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_unif_strains \
    -fiber_cilin 0.4,-1.0,+1.0 \
    -fiber_nx $i \
    -fiber_ny $i \
    -print_vtu > micro_$i.out

    awk -v i_awk=$i '
    /Constitutive/ {
      printf "%d ",i_awk;
      for(j=0;j<6;j++){
	getline;
	for(i=1;i<=6;i++){
	  printf "%e ",$i;
	}
      }
      printf "\n";
    }' micro_$i.out >> run_4_a.dat

    if [ -d "run_4/rve_${i}_a" ]; then
       rm -rf run_4/rve_${i}_a/*
    else
       mkdir run_4/rve_${i}_a
    fi

    mv micro_*  run_4/rve_${i}_a/.
    echo "ok"

done

mv run_4_a.dat run_4/.

#---------------------------------------------------------------------- 

for i in $(seq 1 5);
do

    echo "case $i"
    ./mpirun -np $NM ../../micro/micro \
    -input ex2.spu \
    -mesh meshes/rve_3/rve_$i.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_unif_strains \
    -fiber_cilin 0.4,-0.5,+0.5 \
    -fiber_nx $i \
    -fiber_ny $i \
    -print_vtu > micro_$i.out

    awk -v i_awk=$i '
    /Constitutive/ {
      printf "%d ",i_awk;
      for(j=0;j<6;j++){
	getline;
	for(i=1;i<=6;i++){
	  printf "%e ",$i;
	}
      }
      printf "\n";
    }' micro_$i.out >> run_4_b.dat

    if [ -d "run_4/rve_${i}_b" ]; then
       rm -rf run_4/rve_${i}_b/*
    else
       mkdir run_4/rve_${i}_b
    fi

    mv micro_*  run_4/rve_${i}_b/.
    echo "ok"

done

mv run_4_b.dat run_4/.

#---------------------------------------------------------------------- 

for i in $(seq 1 5);
do

    echo "case $i"
    ./mpirun -np $NM ../../micro/micro \
    -input ex2.spu \
    -mesh meshes/rve_3/rve_$i.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_unif_strains \
    -fiber_cilin 0.4,-0.0,+0.0 \
    -fiber_nx $i \
    -fiber_ny $i \
    -print_vtu > micro_$i.out

    awk -v i_awk=$i '
    /Constitutive/ {
      printf "%d ",i_awk;
      for(j=0;j<6;j++){
	getline;
	for(i=1;i<=6;i++){
	  printf "%e ",$i;
	}
      }
      printf "\n";
    }' micro_$i.out >> run_4_c.dat

    if [ -d "run_4/rve_${i}_c" ]; then
       rm -rf run_4/rve_${i}_c/*
    else
       mkdir run_4/rve_${i}_c
    fi

    mv micro_*  run_4/rve_${i}_c/.
    echo "ok"

done

mv run_4_c.dat run_4/.
