#!/bin/bash
# 强扩展性测试脚本：总问题规模固定为 2^22，进程数变化

for P in 1 2 4 8; do
    echo "Running Strong Scaling Test with P=$P"
    mpirun -np $P ./sum > strong_scaling_P${P}.log
done