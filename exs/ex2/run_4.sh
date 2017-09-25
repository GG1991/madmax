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

size= ( 10 )

for i in ${size[@]};
do

  ./mpirun -np $NM ../../micro/micro \
    -input ex2.spu \
    -mesh meshes/rve_s/rve_$i.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -fiber_cilin 0.4,0.0,0.0 \
    -homo_ld \
    -print_vtu > micro.out

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
    }' micro.out >> run_4.dat

done

