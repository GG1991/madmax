#!/bin/bash
./mpirun -np 1 ../../macro/macro ex1.spu \
         -pc_type lu \
	 -p_vtk   2  \
       : -np 1 ../../micro/micro ex1.spu
