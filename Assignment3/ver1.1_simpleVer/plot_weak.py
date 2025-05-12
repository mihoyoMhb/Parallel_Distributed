import matplotlib.pyplot as plt

# 原始数据（仅处理随机输入部分，即"Sorted"其实为Random）
process_counts = [1, 2, 4, 8, 16]

# 每种策略的运行时间（单位：秒）
strategy_1_times = [21.503977, 23.189609, 26.580503, 31.177770, 32.838663]
strategy_2_times = [21.472940, 23.177135, 26.478694, 31.103781, 32.958514]
strategy_3_times = [21.505544, 23.251905, 26.619421, 31.106057, 32.865768]

# 效率计算（基于每种策略的1进程时间）
efficiency_1 = [strategy_1_times[0] / t for t in strategy_1_times]
efficiency_2 = [strategy_2_times[0] / t for t in strategy_2_times]
efficiency_3 = [strategy_3_times[0] / t for t in strategy_3_times]

# 使用默认的 matplotlib 风格重新绘图
plt.style.use("default")

# 图1：运行时间 vs 进程数（弱缩放）
plt.figure(figsize=(10, 6))
plt.plot(process_counts, strategy_1_times, marker='o', label='Pivot Strategy 1')
plt.plot(process_counts, strategy_2_times, marker='s', label='Pivot Strategy 2')
plt.plot(process_counts, strategy_3_times, marker='^', label='Pivot Strategy 3')
plt.title("Weak Scaling Runtime (Random Input)")
plt.xlabel("Number of Processes")
plt.ylabel("Runtime (seconds)")
plt.xticks(process_counts)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# 图2：效率 vs 进程数
plt.figure(figsize=(10, 6))
plt.plot(process_counts, efficiency_1, marker='o', label='Pivot Strategy 1')
plt.plot(process_counts, efficiency_2, marker='s', label='Pivot Strategy 2')
plt.plot(process_counts, efficiency_3, marker='^', label='Pivot Strategy 3')
plt.title("Weak Scaling Efficiency (Random Input)")
plt.xlabel("Number of Processes")
plt.ylabel("Efficiency")
plt.xticks(process_counts)
plt.ylim(0, 1.1)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# 逆序输入的运行时间（单位：秒）
backward_1_times = [6.366224, 7.675448, 10.374088, 13.713047, 16.784279]
backward_2_times = [6.366550, 7.412714, 9.040388, 11.684343, 13.220409]
backward_3_times = [6.367958, 7.635895, 9.352339, 11.685081, 13.024088]

# 效率计算（基于每种策略的1进程时间）
eff_b1 = [backward_1_times[0] / t for t in backward_1_times]
eff_b2 = [backward_2_times[0] / t for t in backward_2_times]
eff_b3 = [backward_3_times[0] / t for t in backward_3_times]

# 图1：运行时间 vs 进程数（加入逆序输入）
plt.figure(figsize=(10, 6))
plt.plot(process_counts, strategy_1_times, marker='o', label='Random - Pivot 1')
plt.plot(process_counts, strategy_2_times, marker='s', label='Random - Pivot 2')
plt.plot(process_counts, strategy_3_times, marker='^', label='Random - Pivot 3')
plt.plot(process_counts, backward_1_times, marker='o', linestyle='--', label='Backwards - Pivot 1')
plt.plot(process_counts, backward_2_times, marker='s', linestyle='--', label='Backwards - Pivot 2')
plt.plot(process_counts, backward_3_times, marker='^', linestyle='--', label='Backwards - Pivot 3')
plt.title("Weak Scaling Runtime: Random vs Backwards Input")
plt.xlabel("Number of Processes")
plt.ylabel("Runtime (seconds)")
plt.xticks(process_counts)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# 图2：效率 vs 进程数（加入逆序输入）
plt.figure(figsize=(10, 6))
plt.plot(process_counts, efficiency_1, marker='o', label='Random - Pivot 1')
plt.plot(process_counts, efficiency_2, marker='s', label='Random - Pivot 2')
plt.plot(process_counts, efficiency_3, marker='^', label='Random - Pivot 3')
plt.plot(process_counts, eff_b1, marker='o', linestyle='--', label='Backwards - Pivot 1')
plt.plot(process_counts, eff_b2, marker='s', linestyle='--', label='Backwards - Pivot 2')
plt.plot(process_counts, eff_b3, marker='^', linestyle='--', label='Backwards - Pivot 3')
plt.title("Weak Scaling Efficiency: Random vs Backwards Input")
plt.xlabel("Number of Processes")
plt.ylabel("Efficiency")
plt.xticks(process_counts)
plt.ylim(0, 1.1)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
