#!/bin/bash
#SBATCH -A uppmax2025-2-247
#SBATCH -n 4
#SBATCH -c 1
#SBATCH --time=00:10:00
#SBATCH --output=weak_4_%j.out
#SBATCH --error=weak_4_%j.err

module load openmpi/5.0.5

make
mpirun ./sum 2097152
