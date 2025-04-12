#!/bin/bash
#SBATCH -A uppmax2025-2-247
#SBATCH -n 16
#SBATCH -c 1
#SBATCH --time=00:20:00
#SBATCH --output=weak_scaling_iters_%j.out
#SBATCH --error=weak_scaling_iters_%j.err

module load openmpi/5.0.5
make clean
make

# 固定大文件
# INPUTFILE="/proj/uppmax2025-2-247/A2/input1000000.txt"
INPUTFILE="/home/mihoyohb/Datas/input1000000.txt"
OUTPUTFILE="/dev/null"

# 基准迭代次数
BASE_ITER=1000

# 1 rank => ITER = 1 × BASE_ITER
echo "1 rank, iteration=$BASE_ITER"
mpirun -n 1 ./stencil $INPUTFILE $OUTPUTFILE $BASE_ITER

# 2 ranks => ITER = 2 × BASE_ITER
IT2=$((BASE_ITER * 2))
echo "2 ranks, iteration=$IT2"
mpirun -n 2 ./stencil $INPUTFILE $OUTPUTFILE $IT2

# 4 ranks => ITER = 4 × BASE_ITER
IT4=$((BASE_ITER * 4))
echo "4 ranks, iteration=$IT4"
mpirun -n 4 ./stencil $INPUTFILE $OUTPUTFILE $IT4

# 8 ranks => ITER = 8 × BASE_ITER
IT8=$((BASE_ITER * 8))
echo "8 ranks, iteration=$IT8"
mpirun -n 8 ./stencil $INPUTFILE $OUTPUTFILE $IT8

# 16 ranks => ITER = 16 × BASE_ITER
IT16=$((BASE_ITER * 16))
echo "16 ranks, iteration=$IT16"
mpirun -n 16 ./stencil $INPUTFILE $OUTPUTFILE $IT16
