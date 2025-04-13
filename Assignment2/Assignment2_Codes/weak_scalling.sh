#!/bin/bash 
#SBATCH -A uppmax2025-2-247
#SBATCH -n 16
#SBATCH -c 1
#SBATCH --time=01:00:00
#SBATCH --output=weak_scaling_iters_%j.out
#SBATCH --error=weak_scaling_iters_%j.err

module load gcc/14.2.0
module load openmpi/5.0.5
make clean
make

# 基准迭代次数
BASE_ITER=1000

# 输入文件列表
INPUTFILES=(
    "/proj/uppmax2025-2-247/A2/input1000000.txt"
    "/proj/uppmax2025-2-247/A2/input2000000.txt"
)

# 每个配置重复测试的次数
REPEAT=5

# 对每个输入文件进行 weak scaling 测试
for INPUTFILE in "${INPUTFILES[@]}"; do
    FILESIZE=$(basename "$INPUTFILE" | sed 's/input\([0-9]*\)\.txt/\1/')
    echo "========== Testing input size: $FILESIZE =========="
    OUTPUTFILE="/dev/null"

    for RANK in 1 2 4 6 8 10 12 14 16; do
        ITER=$((BASE_ITER * RANK))
        echo "$RANK ranks, iteration=$ITER"

        sum=0
        for ((i=1; i<=REPEAT; i++)); do
            echo "  Run $i:"
            TIME=$(mpirun --bind-to none -n $RANK ./stencil $INPUTFILE $OUTPUTFILE $ITER | grep -Eo '[0-9]+\.[0-9]+')
            echo "    Time: $TIME s"
            sum=$(awk "BEGIN {print $sum + $TIME}")
        done

        avg=$(awk "BEGIN {print $sum / $REPEAT}")
        echo "  Average time for $RANK ranks: $avg s"
    done
done
