#!/bin/bash
#
# Performs an homogenization for different mesh and prints
# the homogenized constitutive tensor
#

mesh_folder="../../meshes/cube_fiber"
NM=4

if [ -e "convergence_const.dat" ]; then
   rm convergence_const.dat
fi

if [ -e "convergence_stress.dat" ]; then
   rm convergence_stress.dat
fi

sizes=(2 3 4 5 10 15 20 25)
for i in ${sizes[@]};
do

  m4 -DNint_m4=$i $mesh_folder/cube_fiber.geo.m4 > cube_fiber.geo
  gmsh -3 cube_fiber.geo
  rm cube_fiber.geo

  ./mpirun -np $NM ../../micro/micro \
    -input ex1.spu \
    -mesh cube_fiber.msh\
    -mesh_gmsh \
    -ksp_type cg \
    -ksp_rtol 1.0e-13 \
    -ksp_atol 1.0e-19 \
    -pc_type bjacobi \
    -options_left 0 \
    -print_vtu \
    -homo_taylor > micro.out

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
    }' micro.out >> convergence_const.dat

    awk -v i_awk=$i '
    /stress_ave/ {
      printf "%d ",i_awk;
      for(i=1;i<=6;i++){
	printf "%e ",$(i+2);
      }
      printf "\n";
      exit;
    }' micro.out >> convergence_stress.dat

done
