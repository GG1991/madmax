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

names=(                   \
"Read_Elems_of_Mesh"      \
"Partition_Mesh"          \
"Calculate_Ghosts_Nodes"  \
"Reenumerates_Nodes"      \
"Read_Coordinates"        \
"Init_Gauss_Points"       \
"Allocate_Mat_and_Vec"    \
"Set_Displ_on_Bou "       \
"Assembly_Jacobian"       \
"Assembly_Residual"       \
"Solve_Linear_System")

nevents=${#names[@]}

if [ -e "names.dat" ]; then
   rm names.dat
fi
for j in ${names[@]}; do
  awk -v pattern=$j '$0 ~ pattern {printf "%s\n", $5}' macro_trace.0 >> names.dat 
done

for i in `seq 1 $n`; do
  file=$( sed -n ${i}p files.dat )
  if [ -e "trace_$((i-1)).dat" ]; then
     rm trace_$((i-1)).dat
  fi
  for j in ${names[@]}; do
     awk -v pattern=$j '$0 ~ pattern {printf "%e\n", $2}' $file >> trace_$((i-1)).dat 
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

scale_x=5;
width=0.1;
sep=0.1;

awk '
/Read_Elems_of_Mesh/       { printf "%-20s\n","ColorAA"}  
/Partition_Mesh/           { printf "%-20s\n","ColorAA"}
/Calculate_Ghosts_Nodes/   { printf "%-20s\n","ColorAA"}
/Reenumerates_Nodes/       { printf "%-20s\n","ColorAA"}
/Read_Coordinates/         { printf "%-20s\n","ColorAA"}
/Init_Gauss_Points/        { printf "%-20s\n","ColorAA"}
/Allocate_Mat_and_Vec/     { printf "%-20s\n","ColorAA"}
/Set_Displ_on_Bou/         { printf "%-20s\n","VIOLETA"}
/Assembly_Jacobian/        { printf "%-20s\n","BLUE"   }
/Assembly_Residual/        { printf "%-20s\n","YELLOW" }
/Solve_Linear_System/      { printf "%-20s\n","RED"    }
' names.dat > colors.dat

paste colors.dat times.dat > times.dat.aux
mv times.dat.aux times.dat

awk -v width_awk=$width -v sep_awk=$sep -v scale_x_awk=$scale_x ' 
BEGIN{i=1}
{
  color[i] = $1;
  for(j=1;j<=NF;j++){
    array[i,j]=$(j+1);
  }
  i++;
}
END{
  nproc = NF
  for(j=1;j<=nproc;j++){
    ymin = 10 -j * (width_awk+sep_awk) - width_awk ;
    ymax = 10 -j * (width_awk+sep_awk);
    for(i=1;i<=NR/2;i++){
      printf "%d\\filldraw[fill=%s, draw=White] (%e, %e) rectangle (%e,%e);\n"  \
,i,color[2*i-1],array[2*i-1,j]*scale_x_awk, ymin, array[2*i,j]*scale_x_awk,ymax
    }
    printf "\n"
  }
}' times.dat >> trace_half.tex

cat trace_head.tex >  trace_final.tex
cat trace_half.tex >> trace_final.tex
cat trace_tail.tex >> trace_final.tex

pdflatex trace_final.tex
