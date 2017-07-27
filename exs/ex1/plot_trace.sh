#!/bin/bash

n=$( ls -l | grep "macro_trace" | wc -l )
ls -l | grep "macro_trace" | awk '{ print $NF }' > files.dat

names=( "Read_Elems_of_Mesh" "Partition_Mesh" "Calculate_Ghosts_Nodes" "Reenumerates_Nodes" )

for i in `seq 1 $n`; do
  file=$( sed -n ${i}p files.dat )
  if [ -e "trace_$((i-1)).dat" ]; then
     rm trace_$((i-1)).dat
  fi
  for j in ${names[@]}; do
     awk -v pattern=$j '$0 ~ pattern {print $2}' $file >> trace_$((i-1)).dat 
     echo $j $file
  done
done
