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

if ! [ -d "mic_run_1" ]; then
  mkdir mic_run_1
fi

./mpirun -np $NM ../../micro/micro \
    -input ex2.spu \
    -mesh meshes/rve_4/rve_3.msh \
    -dim 2 \
    -mesh_gmsh \
    -pc_type lu \
    -options_left 0 \
    -homo_us \
    -print_vtu \
    -options_left 0

mv micro_* mic_run_1/.

#-np $NM xterm -e gdb --args
