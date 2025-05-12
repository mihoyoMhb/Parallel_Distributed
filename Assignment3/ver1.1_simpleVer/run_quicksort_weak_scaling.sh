#!/bin/bash
#SBATCH -A uppmax2025-2-247    # 请替换为您的实际项目账户
#SBATCH -N 1                   # 请求 1 个节点
#SBATCH -n 16                  # 为作业请求最多 16 个 MPI rank (核心/进程槽)
#SBATCH -c 1                   # 每个 MPI rank 使用 1 个核心
#SBATCH --time=08:00:00        # 预计运行时间 (弱缩放，大问题规模可能需要更久)
#SBATCH --output=quicksort_weak_scaling_compact_%j.out
#SBATCH --error=quicksort_weak_scaling_compact_%j.err

# 加载所需模块 (版本号请根据您集群的实际情况确认)
module load gcc/14.2.0     # 或者您集群上可用的 gcc 版本
module load openmpi/5.0.5   # 或者您集群上可用的 OpenMPI 版本

# 清理并编译程序
make clean
make

# 定义程序和通用参数
EXECUTABLE="./quicksort"
OUTPUTFILE="/dev/null" # Output of the sort is discarded
BASE_WORKLOAD_PER_PROCESS="125000000"

# 定义MPI进程数数组
declare -a MPI_PROCESSES=(1 2 4 8 16)

# 定义与每个进程数对应的总问题规模
declare -a TOTAL_PROBLEM_SIZES=(
    "125000000"    # For P=1
    "250000000"    # For P=2
    "500000000"    # For P=4
    "1000000000"   # For P=8
    "2000000000"   # For P=16
)

# 定义与每个进程数对应的输入文件片段
declare -a INPUT_FILE_BASENAMES=(
    "input125000000.txt"
    "input250000000.txt"
    "input500000000.txt"
    "input1000000000.txt"
    "input2000000000.txt"
)

# 定义输入文件类型
declare -a INPUT_TYPES=("sorted" "backwards")
declare -a INPUT_TYPE_DESCRIPTIONS=("Sorted" "Backwards") # Simplified descriptions

# 定义主元选择策略数组
declare -a PIVOT_STRATEGIES=(1 2 3)

INPUT_DIR_BASE="/crex/proj/mixed-precision/nobackup/A3/inputs"

echo "Weak Scaling Test - Parallel Quicksort"
echo "User: mihoyoMhb, Date: 2025-05-12 12:34:16 UTC" # Using provided current date
echo "BaseWorkloadPerProcess: $BASE_WORKLOAD_PER_PROCESS"
echo "Parameters: DataType,TotalN,PivotStrategy,Processes,InputFile" # Header for parameter lines
echo "--- Begin Runs ---"

# 最外层循环：遍历输入文件类型 (sorted, backwards)
for type_idx in "${!INPUT_TYPES[@]}"; do
    INPUT_TYPE=${INPUT_TYPES[$type_idx]}
    INPUT_TYPE_DESC=${INPUT_TYPE_DESCRIPTIONS[$type_idx]}

    # 中层循环：遍历主元选择策略
    for PIVOT_STRATEGY in "${PIVOT_STRATEGIES[@]}"; do

        # 内层循环：遍历MPI进程数 (同时也决定了总问题规模)
        for proc_idx in "${!MPI_PROCESSES[@]}"; do
            processes=${MPI_PROCESSES[$proc_idx]}
            CURRENT_TOTAL_N=${TOTAL_PROBLEM_SIZES[$proc_idx]}
            FILE_BASENAME=${INPUT_FILE_BASENAMES[$proc_idx]}
            
            INPUTFILE=""
            if [ "$INPUT_TYPE" == "sorted" ]; then
                INPUTFILE="$INPUT_DIR_BASE/$FILE_BASENAME"
            else # backwards
                INPUTFILE="$INPUT_DIR_BASE/backwards/input_backwards${FILE_BASENAME#input}"
            fi

            # 打印紧凑的参数行
            echo "RunParams,$INPUT_TYPE_DESC,$CURRENT_TOTAL_N,$PIVOT_STRATEGY,$processes,$INPUTFILE"
            
            # 实际执行命令 (您的程序应将其执行时间打印到标准输出)
            mpiexec -n $processes $EXECUTABLE $INPUTFILE $OUTPUTFILE $PIVOT_STRATEGY
            
        done
    done
done

echo "--- All Runs Finished ---"