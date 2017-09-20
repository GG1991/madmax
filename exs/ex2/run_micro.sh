#!/bin/bash

NM=1

#./mpirun -np $NM ../../micro/micro \
#    -input ex2.spu \
#    -mesh cube_fiber_2d.msh \
#    -dim 2 \
#    -mesh_gmsh \
#    -pc_type lu \
#    -options_left 0 \
#    -homo_unif_strains \
#    -print_vtu

./mpirun -np $NM ../../micro/micro \
    -np $Nm ../../micro/micro \
    -input ex2.spu \
    -mesh cube_2d.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_unif_strains \
    -print_vtu \
    -fiber_cilin 0.4,-0.0,+0.0 \
    -fiber_nx 1 \
    -fiber_ny 1 \
    -options_left 0

#-np $NM xterm -e gdb --args
