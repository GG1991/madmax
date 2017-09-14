#!/bin/bash
#
# Performs an homogenization for different mesh and prints
# the homogenized constitutive tensor
#

NM=1


function cube_fiber_2d {

if [ -e "conv_c.dat" ]; then
   rm conv_c.dat
fi

if [ -e "conv_stress.dat" ]; then
   rm conv_stress.dat
fi
mesh_folder="../../meshes/cube_fiber"
sizes=(7 8 9 10 15 )
for i in ${sizes[@]};
do

  m4 -DNint_m4=$i $mesh_folder/cube_fiber_2d_analysis.geo.m4 > cube_fiber_2d_analysis.geo
  gmsh -2 cube_fiber_2d_analysis.geo
  rm cube_fiber_2d_analysis.geo

  ./mpirun -np $NM ../../micro/micro \
    -input ex2.spu \
    -mesh cube_fiber_2d_analysis.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_ld \
    -print_vtu > micro.out
    #-fiber_cilin 0.4 \

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
    }' micro.out >> conv_c.dat

    awk -v i_awk=$i '
    /stress_ave/ {
      printf "%d ",i_awk;
      for(i=1;i<=6;i++){
	printf "%e ",$(i+2);
      }
      printf "\n";
      exit;
    }' micro.out >> conv_stress.dat

done
}

cube_fiber_2d
