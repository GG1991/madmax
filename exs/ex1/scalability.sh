#!/bin/bash

./set_problem.sh $1
np=( 1 2 4 )
for i in ${np[@]};do

exec_mac="../../macro/macro ex1.spu -log_trace macro_trace"
exec_mic="../../micro/micro ex1.spu"
echo "mpirun -np $i "$exec_mac" : -np 1 "$exec_mic""
eval  mpirun -np $i "$exec_mac" : -np 1 "$exec_mic"

mv macro_trace.0 macro_trace.0_${i}p

done


