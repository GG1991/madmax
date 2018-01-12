#!/bin/bash

MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec" 

solvers=('-ksp_type cg -pc_type jacobi -ksp_atol 1.0e-9 -ksp_rtol 1.0e-9' '-ksp_type cg -ksp_atol 1.0e-9 -ksp_rtol 1.0e-9' '-pc_type lu')
sizes=(50 100 150 200 250 300 350)

for sol in $(seq 0 2); do
for siz in $(seq 0 6); do

$MPIEXEC -np 1  ../../micro/micro \
    -struct_n ${sizes[$siz]},${sizes[$siz]} \
    -dim 2 \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3" \
    -micro_struct "fiber_cilin 3.0 3.0 1 1 0.75 0.0 0.0" \
    ${solvers[$sol]} \
    -homo_us \
    -print_pvtu > "test_"$sol"_"$siz".dat"

echo $sol, $siz, ${solvers[$sol]}

done
done

rm -f sizes.dat
for siz in $(seq 1 ${#sizes[@]}); do
  echo ${sizes[$((siz-1))]}" " >> sizes.dat
done

rm -f times_sol_0.dat times_sol_1.dat times_sol_2.dat

for sol in $(seq 0 2); do
for siz in $(seq 0 6); do

  awk 'BEGIN{i=0}/time sol/{if (i == 1) {print $4; exit; } else { i++; }}' "test_"$sol"_"$siz".dat" >> "times_sol_"$sol".dat"

done
cut -c1-8 "times_sol_"$sol".dat" > "times_sol_"$sol"_cut.dat"
done

paste sizes.dat times_sol_0_cut.dat > times_sol_cg_jp.dat
paste sizes.dat times_sol_1_cut.dat > times_sol_cg.dat
paste sizes.dat times_sol_2_cut.dat > times_sol_lu.dat
