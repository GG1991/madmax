#!/bin/bash

#break_mac=( 'mac_main.c:xx' ) 
#break_mic=( 'mic_main.c:xx' ) 
#break_mac=( 'spu_mesh.c:xx' ) 
#break_mic=( 'spu_mesh.c:xx' ) 
#break_mac=( 'macmic.c:xx' ) 
#break_mic=( 'macmic.c:xx' ) 
#break_mac=( 'mac_alloc.c:xx' ) 
#break_mac=( 'spu_assembly.c:xx' ) 
#break_mac=( 'spu_boundary.c:xx' ) 

NM=1
Nm=1


# BREAKPOINTS
for i in ${break_mac[@]}
do
  exopt_mac="$exopt_mac -ex 'break $i' "
done
exopt_mac+="-ex 'r'"

for i in ${break_mic[@]}
do
  exopt_mic="$exopt_mic -ex 'break $i' "
done
exopt_mic+="-ex 'r'"

gdbcomm_mac="gdb $exopt_mac --args  ../../macro/macro \
	     -input ex1.spu \
	     -mesh ../../meshes/cube_unif/cube.msh \
	     -coupl \
	     -mesh_gmsh \
	     -ksp_type cg \
	     -ksp_rtol 1.0e-13 \
	     -ksp_atol 1.0e-11 \
	     -coupl \
	     -tf 1.0 \
	     -dt 1.0 \
	     -options_left 0 \
	     -testcomm 1"

gdbcomm_mic="gdb $exopt_mic --args ../../micro/micro \
	     -input ex1.spu \
	     -mesh ../../meshes/cube_unif/cube.msh \
	     -mesh_gmsh \
	     -ksp_type cg \
	     -ksp_rtol 1.0e-13 \
	     -ksp_atol 1.0e-11 \
	     -options_left 0 \
	     -homo_exp \
	     -coupl "

./mpirun -np $NM xterm -e "$gdbcomm_mac" : -np $Nm xterm -e "$gdbcomm_mic"
