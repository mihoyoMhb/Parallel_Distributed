#!/bin/bash
#SBATCH -A uppmax2025-2-247
#SBATCH -n 16            # 申请最多 16 个 MPI rank（进程）
#SBATCH -c 1
#SBATCH --time=00:20:00
#SBATCH --output=strong_scaling_%j.out
#SBATCH --error=strong_scaling_%j.err

module load openmpi/5.0.5

make clean
make

# 输入文件路径
# INPUTFILE="/proj/uppmax2025-2-247/A2/input8000000.txt"
INPUTFILE="/home/mihoyohb/Datas/input8000000.txt"
# 输出文件使用 /dev/null，这样程序不会实际写入数据
OUTPUTFILE="/dev/null"
# 迭代次数，根据实际需求进行调整，比如 10 次
ITER=100

echo "Running 1 MPI rank..."
mpirun -n 1 ./stencil $INPUTFILE $OUTPUTFILE $ITER

echo "Running 2 MPI ranks..."
mpirun -n 2 ./stencil $INPUTFILE $OUTPUTFILE $ITER

echo "Running 4 MPI ranks..."
mpirun -n 4 ./stencil $INPUTFILE $OUTPUTFILE $ITER

echo "Running 8 MPI ranks..."
mpirun -n 8 ./stencil $INPUTFILE $OUTPUTFILE $ITER

echo "Running 16 MPI ranks..."
mpirun -n 16 ./stencil $INPUTFILE $OUTPUTFILE $ITER
