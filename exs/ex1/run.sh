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

./mpirun -np 2 ../../macro/macro ex1.spu \
         -ksp_type cg  \
	 -p_vtk   2  \
       : -np 1 ../../micro/micro ex1.spu
