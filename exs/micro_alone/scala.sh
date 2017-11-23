#!/bin/bash

if [ -e "times.dat" ]
then
  rm times.dat
fi

echo "running"
for i in $(seq 1 10) 
do
  echo "$i proc"
  ./run.sh $i > out.dat
  mv micro_trace.dat micro_trace_$i.dat
  echo "$i $(sed '$!d' micro_trace_$i.dat | awk '{print $2}')" >> times.dat

done
