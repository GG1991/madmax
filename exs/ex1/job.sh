#!/bin/bash
#SBATCH --job-name="sputnik_g"
#SBATCH --time=00-00:02:00
#SBATCH --output=sputnik_%j.out
#SBATCH --error=sputnik_%j.out
#SBATCH --ntasks=1025
mpirun -np 1024 ../../macro/macro ex1.spu : -np 1 ../../micro/micro ex1.spu
