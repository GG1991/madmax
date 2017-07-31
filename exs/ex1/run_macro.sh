#!/bin/bash

NM=2

./mpirun -np $NM ../../macro/macro ex1.spu \
    -ksp_type cg                    \
    -ksp_rtol 1.0e-13
#    -coupl 1                        \
