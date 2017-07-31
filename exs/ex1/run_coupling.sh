#!/bin/bash
#./mpirun -np 1 ../../macro/macro ex1.spu \
#         -ksp_type cg  \
#	 -p_vtk   2  \
#       : -np 1 ../../micro/micro ex1.spu

#./mpirun -np 1 ../../macro/macro ex1.spu \
#         -ksp_type cg  \
#         -pc_type lu  \
#	 -p_vtk   2  \
#       : -np 1 ../../micro/micro ex1.spu
NM=2
Nm=2

./mpirun -np $NM ../../macro/macro ex1.spu      \
         -ksp_type cg                           \
	 -ksp_rtol 1.0e-13                      \
	 -log_trace macro_trace                 \
	 -p_vtk   2                             \
	 -coupl                                 \
	 -options_left 0                        \
	 -mesh ../../meshes/cube_unif/cube.msh  \
       : -np $Nm ../../micro/micro ex1.spu      \
	 -coupl                                 \
	 -options_left 0                        \
	 -mesh ../../meshes/cube_unif/cube.msh

#      -ksp_atol <abstol> 	- Sets abstol
#      -ksp_rtol <rtol> 	- Sets rtol
#      -ksp_divtol <dtol> 	- Sets dtol
#      -ksp_max_it <maxits> 	- Sets maxits 
