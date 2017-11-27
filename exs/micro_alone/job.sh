#!/bin/bash
# the script is design to run on the following folder and for being launched with sbatch ../job.sh

#SBATCH --job-name="micro_1"
#SBATCH --time=00-00:01:30
#SBATCH --output=micro_%j.out
#SBATCH --error=micro_%j.out
#SBATCH --nodes=1

MPIEXEC="mpiexec" 
EXECUTABLE="../../../micro/micro"

srun $EXECUTABLE \
    -struct_n      75,75 \
    -dim           2       \
    -material      "MATRIX TYPE_0 1.0e7 1.0e6 0.3","FIBER TYPE_0 1.0e7 1.0e7 0.3" \
    -micro_struct  "fiber_line 3.0 3.0 2 9 9 0.785398 2.35619 0.4 0.4 0.2 0.2 0.0 0.0" \
    -pc_type       jacobi  \
    -ksp_type      cg      \
    -print_vtu             \
    -homo_us               \
    -options_left  0
