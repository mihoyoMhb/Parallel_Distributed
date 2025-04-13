import matplotlib.pyplot as plt


# 定义每个输入规模对应的运行时间数据（均为迭代 2500 次）
scaling_data = {
    "1M": {
        "ranks": [1, 2, 4, 8, 16],
        "times": [4.136495, 2.127032, 1.226924, 0.780439, 0.754196]
    },
    "2M": {
        "ranks": [1, 2, 4, 8, 16],
        "times": [8.250402, 4.197815, 2.341748, 1.426417, 1.188866]
    },
    "4M": {
        "ranks": [1, 2, 4, 8, 16],
        "times": [16.558072, 8.450844, 4.467812, 2.634149, 2.206028]
    },
    "8M": {
        "ranks": [1, 2, 4, 8, 16],
        "times": [29.721115, 14.803920, 8.112024, 5.294519, 4.523037]
    }
}

plt.figure(figsize=(10, 6))
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
linestyles = ['-', '--', '-.', ':']

# 画加速比曲线
for i, (label, d) in enumerate(scaling_data.items()):
    base_time = d["times"][0]
    speedups = [base_time / t for t in d["times"]]
    plt.plot(d["ranks"], speedups, marker='o', label=f'{label} points', color=colors[i], linestyle=linestyles[i])

# 理想加速比线
ranks = scaling_data["1M"]["ranks"]
plt.plot(ranks, ranks, 'k--', label='Ideal speedup')

plt.xlabel("Number of MPI ranks")
plt.ylabel("Speedup")
plt.title("Strong Scaling with 2500 Iterations (Different Input Sizes)")
plt.legend()
plt.grid(True)
plt.xticks(ranks)
plt.tight_layout()
plt.show()


# Ranks and average times
ranks = [1, 2, 4, 6, 8, 10, 12, 14, 16]
time_1m = [1.12886, 1.13886, 1.21228, 1.29563, 1.30558, 1.34201, 1.37492, 1.39909, 1.46107]
time_2m = [2.62112, 2.24634, 2.26888, 2.39058, 2.46892, 2.5599, 2.63926, 2.70796, 2.79568]

# Efficiency = T1 / Tp
eff_1m = [time_1m[0] / t for t in time_1m]
eff_2m = [time_2m[0] / t for t in time_2m]

# Plot combined time figure
plt.figure(figsize=(9, 5))
plt.plot(ranks, time_1m, marker='o', label='1M Elements')
plt.plot(ranks, time_2m, marker='s', label='2M Elements')
plt.title('Weak Scaling: Average Computation Time')
plt.xlabel('Number of MPI Ranks')
plt.ylabel('Average Time (s)')
plt.grid(True)
plt.xticks(ranks)
plt.legend()
plt.tight_layout()
plt.show()

# Plot combined efficiency figure
plt.figure(figsize=(9, 5))
plt.plot(ranks, eff_1m, marker='o', label='1M Elements')
plt.plot(ranks, eff_2m, marker='s', label='2M Elements')
plt.title('Weak Scaling: Efficiency (T1 / Tp)')
plt.xlabel('Number of MPI Ranks')
plt.ylabel('Efficiency')
plt.grid(True)
plt.xticks(ranks)
plt.legend()
plt.tight_layout()
plt.show()
