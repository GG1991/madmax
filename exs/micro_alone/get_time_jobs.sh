#!/bin/bash

if [ -e "times.dat" ]
then
  rm times.dat
fi

for i in $(seq 1 24) 
do
  if [ ! -d "micro_$i" ]
  then
    echo "directory micro_$i not found ... bye"
    exit 1
  fi

  cd micro_$i
  echo "$i $(sed '$!d' micro_trace.dat | awk '{print $2}')" >> ../times.dat
  cd ..

done

