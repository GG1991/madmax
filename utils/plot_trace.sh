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

width=0.1;
sep=0.1;
for ip in `seq 1 $n`; do
 y_min=`echo "10-${ip}*${width}-${sep}"|bc -l`
 y_max=`echo "10-${ip}*${width}-${sep}+${width}"|bc -l`
 for j in `seq 1 $nevents`; do
 
 k=`echo "2*${j}-1"|bc -l`
 eval sed -n ${k}p times.dat > line.dat
 start_num=$(awk -v pattern=${ip} '{print $pattern}' line.dat )

 k=`echo "2*${j}"|bc -l`
 eval sed -n ${k}p times.dat > line.dat
 end_num=$(awk '{print $1}' line.dat )

 echo $start_num $end_num
 
 echo "\filldraw[fill=ColorAA, draw=White] \
 (${start_num}, ${y_min}) rectangle (${end_num}, ${y_max});">> trace_half.tex
 
 done
done


cat trace_head.tex >  trace_final.tex
cat trace_half.tex >> trace_final.tex
cat trace_tail.tex >> trace_final.tex

#pdflatex trace_final.tex
