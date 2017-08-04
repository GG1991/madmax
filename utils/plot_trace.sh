#!/bin/bash
#
# This script takes PETSc traces generated with 
# -log_trace option and draws a tikz picture
#
#  Author : Guido Giuntoli
#  Start date: 27-07-2017
#


n=$( ls -l | grep "macro_trace" | wc -l )
ls -l | grep "macro_trace" | awk '{ print $NF }' > files.dat

names=(                  \
"Read_Elems_of_Mesh"     \
"Partition_Mesh"         \
"Calculate_Ghosts_Nodes" \
"Reenumerates_Nodes" )
nevents=${#names[@]}

if [ -e "names.dat" ]; then
   rm names.dat
fi
for j in ${names[@]}; do
  echo $j >> names.dat 
done

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

cat trace_0.dat > aux1
for i in `seq 2 $n`; do
     paste aux1 trace_$((i-1)).dat > aux2
     cat aux2 > aux1
done
mv aux2 times.dat	
rm aux1

# Armamos los rectangulitos a partir de cada fila de times.dat
if [ -e "trace_half.tex" ]; then
   rm trace_half.tex
fi

echo "\filldraw[fill=ColorAA, draw=White] (0.0, 0.0) rectangle (1.0, 1.0);"> trace_half.tex

cat trace_head.tex >  trace_final.tex
cat trace_half.tex >> trace_final.tex
cat trace_tail.tex >> trace_final.tex

pdflatex trace_final.tex
