#!/bin/bash

#break_mac=( 'mac_main.c:38' ) 
#break_mic=( 'mic_main.c:32' ) 
#break_mac=( 'mac_comm.c:100' ) 
#break_mic=( 'mic_comm.c:104' ) 
break_mac=( 'mac_comm.c:107' ) 
break_mic=( 'mic_comm.c:111' ) 


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

gdbcomm_mac="gdb $exopt_mac --args  ../../macro/macro ../../macro/exs/ex1/ex1.mac"
gdbcomm_mic="gdb $exopt_mic --args  ../../micro/micro ../../macro/exs/ex1/ex1.mac"

mpirun -np 1 xterm -e "$gdbcomm_mac" : -np 1 xterm -e "$gdbcomm_mic"
