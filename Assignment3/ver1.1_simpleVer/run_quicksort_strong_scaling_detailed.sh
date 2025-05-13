#!/bin/bash
#SBATCH -A uppmax2025-2-247  # 请替换为您的实际项目账户
#SBATCH -N 1                 # 请求 1 个节点
#SBATCH -n 16                # 为作业请求最多 16 个 MPI rank (进程槽)
#SBATCH -c 1                 # 每个 MPI rank 使用 1 个核心
#SBATCH --time=05:00:00      # 预计运行时间 (由于增加了测试文件，时间已增加)
#SBATCH --output=quicksort_strong_scaling_multi_size_%j.out
#SBATCH --error=quicksort_strong_scaling_multi_size_%j.err

# 加载所需模块 (版本号请根据您集群的实际情况确认)
module load gcc/14.2.0   # 或者您集群上可用的 gcc 版本
module load openmpi/5.0.5 # 或者您集群上可用的 OpenMPI 版本

# 清理并编译程序
make clean
make

# 定义程序和通用参数
EXECUTABLE="./quicksort" # 脚本中已有的可执行文件名
OUTPUTFILE="/dev/null"

# 定义输入文件数组
declare -a INPUT_FILES=(
    "/crex/proj/mixed-precision/nobackup/A3/inputs/input250000000.txt"
    "/crex/proj/mixed-precision/nobackup/A3/inputs/backwards/input_backwards250000000.txt"
    "/crex/proj/mixed-precision/nobackup/A3/inputs/input125000000.txt"
    "/crex/proj/mixed-precision/nobackup/A3/inputs/backwards/input_backwards125000000.txt"
)
declare -a INPUT_FILE_DESCRIPTIONS=( # 用于日志输出
    "Random_Input_N250M"
    "Backwards_Input_N250M"
    "Random_Input_N125M"
    "Backwards_Input_N125M"
)
declare -a PROBLEM_SIZES=( # 对应每个输入文件的问题规模
    "250000000"
    "250000000"
    "125000000"
    "125000000"
)

# 定义主元选择策略数组
declare -a PIVOT_STRATEGIES=(1 2 3)

# 定义MPI进程数数组
declare -a MPI_PROCESSES=(1 2 4 8 16)


echo "======================================================================"
echo "Parallel Quicksort - Detailed Strong Scaling Test (Multiple Sizes)"
echo "Output File: $OUTPUTFILE"
echo "Current User: mihoyoMhb"
echo "Submission Date (UTC): 2025-05-12 11:43:40" # User provided time
echo "Actual Run Date: $(date)"
echo "======================================================================"
echo ""

# 外层循环：遍历输入文件
for i in "${!INPUT_FILES[@]}"; do
    INPUTFILE=${INPUT_FILES[$i]}
    INPUT_DESC=${INPUT_FILE_DESCRIPTIONS[$i]}
    CURRENT_PROBLEM_SIZE=${PROBLEM_SIZES[$i]}

    echo "######################################################################"
    echo "Testing Input File: $INPUTFILE"
    echo "Description: $INPUT_DESC"
    echo "Problem Size (N): $CURRENT_PROBLEM_SIZE elements"
    echo "######################################################################"
    echo ""

    # 中层循环：遍历主元选择策略
    for PIVOT_STRATEGY in "${PIVOT_STRATEGIES[@]}"; do
        echo "    =============================================================="
        echo "    Input: $INPUT_DESC, Pivot Strategy: $PIVOT_STRATEGY"
        echo "    =============================================================="
        echo ""

        # 内层循环：遍历MPI进程数
        for processes in "${MPI_PROCESSES[@]}"; do
            echo "        ------------------------------------------------------"
            echo "        Input: $INPUT_DESC, Problem Size: $CURRENT_PROBLEM_SIZE, Pivot: $PIVOT_STRATEGY, MPI Ranks: $processes"
            echo "        Command: mpiexec -n $processes $EXECUTABLE $INPUTFILE $OUTPUTFILE $PIVOT_STRATEGY"
            
            # 实际执行命令，您的程序应该会将运行时间打印到标准输出
            mpiexec -n $processes $EXECUTABLE $INPUTFILE $OUTPUTFILE $PIVOT_STRATEGY
            
            echo "        Finished run."
            echo "        ------------------------------------------------------"
            echo ""
        done
        echo "" # Add a blank line after all process counts for a pivot strategy
    done
    echo "" # Add a blank line after all pivot strategies for an input file
done

echo "======================================================================"
echo "All detailed strong scaling tests finished."
echo "======================================================================"