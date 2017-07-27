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

./mpirun -np 4 ../../macro/macro ex1.spu \
         -ksp_type cg                    \
	 -ksp_rtol 1.0e-13               \
	 -log_trace macro_trace          \
	 -p_vtk   2                      \
       : -np 1 ../../micro/micro ex1.spu

#      -ksp_atol <abstol> 	- Sets abstol
#      -ksp_rtol <rtol> 	- Sets rtol
#      -ksp_divtol <dtol> 	- Sets dtol
#      -ksp_max_it <maxits> 	- Sets maxits 
