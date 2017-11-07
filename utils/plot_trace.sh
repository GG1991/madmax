#!/bin/bash
#
#  This script takes the trace generate by "trace.h"
#  and draws a tikz picture
#
#  Author : Guido Giuntoli
#  Start date: 27-07-2017
#

if [ $# -ne 1 ]
then
   echo "Usage ./plot_trace.sh [trace_name]" 1>&2
   exit 1
fi

file=$1

if [ -e "trace_half.tex" ]; then
   rm trace_half.tex
fi

scale_x=5;
width=0.1;
sep=0.1;

awk -v width_awk=$width     \
    -v sep_awk=$sep         \
    -v scale_x_awk=$scale_x \
' 
/ass/{
  for( j = 0 ; j < NF-1 ; j++ )
    start_time[j] = $(j+2)
  getline;
  for( j = 0 ; j < NF-1 ; j++ )
    end_time[j]   = $(j+2)
  for( j = 0 ; j < NF-1 ; j++ ){
      ymin = 10 - j * (width_awk + sep_awk) - width_awk ;
      ymax = 10 - j * (width_awk + sep_awk);
      printf "\\filldraw[fill=%s, draw=White] (%e, %e) rectangle (%e,%e);\n"  \
,"BLUE",start_time[j]*scale_x_awk, ymin, end_time[j]*scale_x_awk,ymax
  }
  printf "\n"
}
/sol/{
  for( j = 0 ; j < NF-1 ; j++ )
    start_time[j] = $(j+2)
  getline;
  for( j = 0 ; j < NF-1 ; j++ )
    end_time[j]   = $(j+2)
  for( j = 0 ; j < NF-1 ; j++ ){
      ymin = 10 - j * (width_awk + sep_awk) - width_awk ;
      ymax = 10 - j * (width_awk + sep_awk);
      printf "\\filldraw[fill=%s, draw=White] (%e, %e) rectangle (%e,%e);\n"  \
,"RED",start_time[j]*scale_x_awk, ymin, end_time[j]*scale_x_awk,ymax
  }
  printf "\n"
}
' $file >> trace_half.tex

cat trace_head.tex >  trace_final.tex
cat trace_half.tex >> trace_final.tex
cat trace_tail.tex >> trace_final.tex

#pdflatex trace_final.tex
