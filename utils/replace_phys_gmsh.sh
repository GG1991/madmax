#!/bin/bash
#
# this script replace "MATRIX_xx" and "FIBER_xx" for "MATRIX" and "FIBER" in a gmsh file
#

if [ $# -ne 1 ];
then
   echo "Usage ./replace_gmsh_phy.sh <gmsh_file>"
   exit 1
fi

ending="_rep"
final_file="${1:0:$(( ${#1} - 4 ))}$ending.msh"

sed -e 's/MATRIX_[0-9]*/MATRIX/g' -e ' s/FIBER_[0-9]*/FIBER/g' $1 > $final_file
echo "generating $final_file"
