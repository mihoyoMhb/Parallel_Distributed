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
