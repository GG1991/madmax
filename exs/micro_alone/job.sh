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
    -material      "MATRIX TYPE_0 1.0e7 1.0e6 0.3","FIBER TYPE_0 1.0e7 1.0e7 0.3" \
    -fiber_cilin   0.03,0.0,0.0,0.0   \
    -dim 2                 \
    -struct_n      150,150 \
    -struct_l      0.1,0.1 \
    -pc_type       jacobi  \
    -ksp_type      cg      \
    -print_vtu             \
    -homo_us               \
    -options_left 0
