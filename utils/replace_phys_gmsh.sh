#!/bin/bash
sed -e 's/MATRIX_[0-9]*/MATRIX/g' -e ' s/FIBER_[0-9]*/FIBER/g' $1 > "${1:0:$(( ${#1} - 4 ))}_a.msh"
echo "generating " "${1:0:$(( ${#1} - 4 ))}_a.msh"

