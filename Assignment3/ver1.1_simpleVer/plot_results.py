import matplotlib.pyplot as plt
import numpy as np

# 实际运行时间数据
ranks = np.array([1, 2, 4, 8, 16])
ideal_speedup = ranks  # 理想加速比

# 按策略和输入组织的数据：{(input_type, strategy): [times]}
timing_data = {
    ("Sorted", 1): [45.572238, 23.730301, 13.187174, 7.463538, 3.821933],
    ("Sorted", 2): [45.579366, 24.163781, 13.225697, 7.457203, 3.812494],
    ("Sorted", 3): [45.579008, 23.585555, 13.179469, 7.466107, 3.814368],
    ("Backwards", 1): [14.404996, 8.263521, 5.398757, 3.394884, 1.981604],
    ("Backwards", 2): [14.380649, 7.940578, 4.706159, 2.832033, 1.538736],
    ("Backwards", 3): [14.383742, 8.230185, 4.832720, 2.892494, 1.549964],
}

# 可视化
plt.figure(figsize=(10, 6))

# 为每组数据绘图
for (input_type, strategy), times in timing_data.items():
    baseline = times[0]
    speedup = [baseline / t for t in times]
    plt.plot(ranks, speedup, marker='o', label=f"{input_type} - Pivot {strategy}")

# 理想加速线
plt.plot(ranks, ideal_speedup, linestyle='--', color='black', label="Ideal Speedup")

# 图形设置
plt.xscale("log", base=2)
plt.yscale("log", base=2)
plt.xticks(ranks, labels=ranks)
plt.yticks(ranks, labels=ranks)
plt.xlabel("Number of MPI Ranks (log2 scale)")
plt.ylabel("Speedup (log2 scale)")
plt.title("Parallel Quicksort - Strong Scaling Speedup (Log-Log Plot)")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.tight_layout()
plt.show()
