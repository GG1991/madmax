#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 
EXECUTABLE="/home/guido/codes/sputnik/micro/micro"

if [ -d c_tang_ave ]; then
  rm c_tang_ave/*;
else
  mkdir c_tang_ave;
fi

#$MPIEXEC -np 1 xterm -e gdb -x file.gdb --args $EXECUTABLE \
$MPIEXEC -np 1 $EXECUTABLE \
    -struct_n 75,75 \
    -dim 2 \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3" \
    -micro_struct "fiber_line 3.0 3.0 1 1 0.0 0.4 1.0 0.0" \
    -pc_type lu \
    -homo_us > c_tang_ave/homo_us.dat;

$MPIEXEC -np 1 $EXECUTABLE \
    -struct_n 75,75 \
    -dim 2 \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3" \
    -micro_struct "fiber_line 3.0 3.0 1 1 0.0 0.4 1.0 0.0" \
    -pc_type lu \
    -homo_tp > c_tang_ave/homo_tp.dat

$MPIEXEC -np 1 $EXECUTABLE \
    -struct_n 75,75 \
    -dim 2 \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3" \
    -micro_struct "fiber_line 3.0 3.0 1 1 0.0 0.4 1.0 0.0" \
    -pc_type lu \
    -homo_ts > c_tang_ave/homo_ts.dat
