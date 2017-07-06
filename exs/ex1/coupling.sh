#!/bin/bash

#break_mac=( 'mac_main.c:66' ) 
#break_mic=( 'mic_main.c:56' ) 
break_mac=( 'spu_mesh.c:136' ) 
break_mic=( 'spu_mesh.c:136' ) 
#break_mac=( 'mac_comm.c:101' ) 
#break_mic=( 'mic_comm.c:105' ) 

NM=10
Nm=1


# BREAKPOINTS
for i in ${break_mac[@]}
do
  exopt_mac="$exopt_mac -ex 'break $i' "
done
exopt_mac+="-ex 'r'"

for i in ${break_mic[@]}
do
  exopt_mic="$exopt_mic -ex 'break $i' "
done
exopt_mic+="-ex 'r'"


if [ "$#" -eq 1 ];then
  if [ "$1" -eq 1 ];then  
   exec_mac="../../macro/macro ex1.spu"
   exec_mic="../../micro/micro ex1.spu"
   echo "mpirun -np $NM "$exec_mac" : -np $Nm "$exec_mic""
   eval  mpirun -np $NM "$exec_mac" : -np $Nm "$exec_mic"
  elif [ "$1" -eq 2 ];then
   exec_val2_mac="valgrind --log-file=\"valgrind_M.out\"  ../../macro/macro ex1.spu"
   exec_val2_mic="valgrind --log-file=\"valgrind_m.out\"  ../../micro/micro ex1.spu"
   echo "mpirun -np $NM "$exec_val2_mac" : -np $Nm "$exec_val2_mic""
   eval  mpirun -np $NM "$exec_val2_mac" : -np $Nm "$exec_val2_mic" 
  elif [ "$1" -eq 3 ];then
   exec_val3_mac="valgrind --leak-check=full ../../macro/macro ex1.spu"
   exec_val3_mic="valgrind --leak-check=full ../../micro/micro ex1.spu"
   echo "mpirun -np 2 "$exec_val3_mac" : -np 2 "$exec_val3_mic""
   eval  mpirun -np 2 "$exec_val3_mac" : -np 2 "$exec_val3_mic" > valgrind3-1.out 2>&1
  elif [ "$1" -eq 4 ];then
   exec_val4_mac="valgrind --leak-check=full ../../macro/macro ex1.spu"
   exec_val4_mic="valgrind --leak-check=full ../../micro/micro ex1.spu"
   echo "mpirun -np 2 "$exec_val4_mac" : -np 2 "$exec_val4_mic""
   eval  mpirun -np 2 "$exec_val4_mac" : -np 2 "$exec_val4_mic"
  fi
else
   gdbcomm_mac="gdb $exopt_mac --args  ../../macro/macro ex1.spu"
   gdbcomm_mic="gdb $exopt_mic --args  ../../micro/micro ex1.spu"
   mpirun -np $NM xterm -e "$gdbcomm_mac" : -np $Nm xterm -e "$gdbcomm_mic"
fi
