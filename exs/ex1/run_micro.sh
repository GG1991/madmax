#!/bin/bash

NM=2

./mpirun -np $NM ../../micro/micro ex1.spu \
    -ksp_type cg                    \
    -ksp_rtol 1.0e-13               \
    -options_left 0
#    -coupl 1
