#!/bin/bash

break_mac=( 'mac_main.c:56' ) 
break_mic=( 'mic_main.c:56' ) 
#break_mac=( 'mac_comm.c:100' ) 
#break_mic=( 'mic_comm.c:104' ) 
#break_mac=( 'mac_comm.c:101' ) 
#break_mic=( 'mic_comm.c:105' ) 


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

gdbcomm_mac="gdb $exopt_mac --args  ../../macro/macro ex1.mac"
gdbcomm_mic="gdb $exopt_mic --args  ../../micro/micro ex1.mac"
exec_mac="../../macro/macro ex1.mac"
exec_mic="../../micro/micro ex1.mac"

#mpirun -np 1 xterm -e "$gdbcomm_mac" : -np 2 xterm -e "$gdbcomm_mic"
if [ "$#" -gt 0 ];then
   echo "mpirun -np 2 "$exec_mac" : -np 2 "$exec_mic""
   eval mpirun -np 2 "$exec_mac" : -np 2 "$exec_mic"
else
   mpirun -np 2 xterm -e "$gdbcomm_mac" : -np 2 xterm -e "$gdbcomm_mic"
fi
