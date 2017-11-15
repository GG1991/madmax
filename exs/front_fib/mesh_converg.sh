#!/bin/bash
#
# Performs an homogenization for different mesh and prints
# the homogenized constitutive tensor
#

NM=1

if [ -e "conv_c.dat" ]; then
   rm conv_c.dat
fi

if [ -e "conv_stress.dat" ]; then
   rm conv_stress.dat
fi

function cube_2d {
mesh_folder="../../meshes/cube_unif"
sizes=(3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 60 65 70 75 80)
for i in ${sizes[@]};
do

  m4 -DN_m4=$i -Dlc_m4=1.5 $mesh_folder/cube_2d_analysis.geo.m4 > cube_2d_analysis.geo
  gmsh -2 cube_2d_analysis.geo
  rm cube_2d_analysis.geo

  ./mpirun -np $NM ../../micro/micro \
    -input ex1_2d.spu \
    -mesh cube_2d_analysis.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_ld \
    -fiber_cilin 0.4 \
    -print_vtu > micro.out

    awk -v i_awk=$i '
    /Constitutive/ {
      printf "%d ",i_awk;
      for(j=0;j<3;j++){
	getline;
	for(i=1;i<=3;i++){
	  printf "%e ",$i;
	}
      }
      printf "\n";
    }' micro.out >> convergence_const_cube_2d.dat

    awk -v i_awk=$i '
    /stress_ave/ {
      printf "%d ",i_awk;
      for(i=1;i<=3;i++){
	printf "%e ",$(i+2);
      }
      printf "\n";
      exit;
    }' micro.out >> convergence_stress_cube_2d.dat

done

}

function cube_fiber_2d {
mesh_folder="../../meshes/cube_fiber"
sizes=(7 8 9 10 15 )
for i in ${sizes[@]};
do

  m4 -DNint_m4=$i $mesh_folder/cube_fiber_2d_analysis.geo.m4 > cube_fiber_2d_analysis.geo
  gmsh -2 cube_fiber_2d_analysis.geo
  rm cube_fiber_2d_analysis.geo

  ./mpirun -np $NM ../../micro/micro \
    -input ex1_2d.spu \
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
}

function rve_size_fiber_2d {
NM=1

for i in $(seq 1 8); do

./mpirun -np $NM ../../micro/micro \
    -input ex1_2d.spu \
    -mesh "../../meshes/cube_fiber/struct_fiber_2d_"$i"_"$i".msh" \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -log_trace micro_trace \
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

}

#
# Creates an NX x NX struct mesh 
# simulates an array of fibers
# is for analize the dependence of the 
# fiber position inside the RVE 
# using LD BC
#
function struc_array {

mesh_folder="../../meshes/cube_unif"
LX=1.5
NX=30
for i in $(seq 1 3);
do

  lx=$(echo "$LX * $i" | bc)
  nx=$(echo "$NX * $i" | bc)
  
  m4 -DN_m4=$nx -Dlc_m4=$lx $mesh_folder/cube_2d_analysis.geo.m4 > cube_2d_analysis.geo
  gmsh -2 cube_2d_analysis.geo
  rm cube_2d_analysis.geo

  ./mpirun -np $NM ../../micro/micro \
    -input ex1_2d.spu \
    -mesh cube_2d_analysis.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_ld \
    -fiber_cilin 0.4,-0.4,+0.4 \
    -fiber_nx $i \
    -fiber_ny $i \
    -print_vtu > micro.out

    awk -v i_awk=$i '
    /Constitutive/ {
      printf "%d ",i_awk;
      for(j=0;j<3;j++){
	getline;
	for(i=1;i<=3;i++){
	  printf "%e ",$i;
	}
      }
      printf "\n";
    }' micro.out >> conv_c.dat

    awk -v i_awk=$i '
    /stress_ave/ {
      printf "%d ",i_awk;
      for(i=1;i<=3;i++){
	printf "%e ",$(i+2);
      }
      printf "\n";
      exit;
    }' micro.out >> conv_stress.dat

done

}
#cube_2d
#cube_fiber_2d
#rve_size_fiber_2d
#struc_array
