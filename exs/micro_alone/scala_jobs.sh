#!/bin/bash

echo "launching jobs"
for i in $(seq 1 24) 
do
  if [ -d "micro_$i" ]
  then
    rm -rf micro_$i/*
  else
    mkdir  micro_$i
  fi

  echo "launching proc $i"
  cd micro_$i
  sbatch --ntasks=$i ../job.sh
  cd ..

done
