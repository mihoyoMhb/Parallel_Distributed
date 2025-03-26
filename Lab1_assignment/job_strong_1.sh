#!/bin/bash
#SBATCH -A uppmax2025-2-247
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=00:10:00
#SBATCH --output=strong_1_%j.out
#SBATCH --error=strong_1_%j.err

module load openmpi/5.0.5

make         # 编译程序（用你的 Makefile）
mpirun ./sum 67108864
