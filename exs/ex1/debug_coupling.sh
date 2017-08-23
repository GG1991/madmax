#!/bin/bash

break_mac=( 'mac_main.c:52' ) 
break_mic=( 'mic_main.c:52' ) 
#break_mac=( 'spu_mesh.c:136' ) 
#break_mic=( 'spu_mesh.c:136' ) 
#break_mac=( 'macmic.c:143' ) 
#break_mic=( 'macmic.c:149' ) 
#break_mac=( 'mac_alloc.c:69' ) 
#break_mac=( 'spu_assembly.c:123' ) 
#break_mac=( 'spu_boundary.c:148' ) 

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
	     -coupl 1 \
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
	     -coupl 1"

./mpirun -np $NM xterm -e "$gdbcomm_mac" : -np $Nm xterm -e "$gdbcomm_mic"
