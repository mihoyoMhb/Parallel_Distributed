import matplotlib.pyplot as plt
import numpy as np

# 数据结构：{input_label: {pivot: [(num_procs, time), ...]}}
data = {
    "Random_Input_N250M": {
        1: [(1, 44.565194), (2, 23.198188), (4, 12.999978), (8, 7.336438), (16, 3.745419)],
        2: [(1, 44.645909), (2, 23.092092), (4, 12.964096), (8, 7.358125), (16, 3.793008)],
        3: [(1, 44.625652), (2, 23.260332), (4, 13.051909), (8, 7.339909), (16, 3.758113)],
    },
    "Backwards_Input_N250M": {
        1: [(1, 13.413998), (2, 7.737786), (4, 5.133247), (8, 3.278565), (16, 1.937009)],
        2: [(1, 13.423428), (2, 7.445135), (4, 4.486369), (8, 2.739015), (16, 1.480628)],
        3: [(1, 13.395346), (2, 7.676731), (4, 4.569595), (8, 2.754746), (16, 1.471175)],
    },
    "Random_Input_N125M": {
        1: [(1, 21.557437), (2, 11.175900), (4, 6.279124), (8, 3.585262), (16, 1.792787)],
        2: [(1, 21.506649), (2, 11.139275), (4, 6.290118), (8, 3.575003), (16, 1.788403)],
        3: [(1, 21.503863), (2, 11.115026), (4, 6.306962), (8, 3.581722), (16, 1.758522)],
    },
    "Backwards_Input_N125M": {
        1: [(1, 6.393453), (2, 3.715719), (4, 2.484808), (8, 1.588038), (16, 0.924150)],
        2: [(1, 6.380606), (2, 3.553009), (4, 2.146394), (8, 1.328349), (16, 0.696411)],
        3: [(1, 6.416172), (2, 3.688495), (4, 2.198494), (8, 1.342028), (16, 0.637620)],
    }
}

# 为每个输入生成图像
for input_name, pivot_data in data.items():
    plt.figure(figsize=(10, 6))
    ideal_line_drawn = False

    for pivot, runs in pivot_data.items():
        procs = np.array([p for p, t in runs])
        times = np.array([t for p, t in runs])
        base_time = times[0]
        speedup = base_time / times

        plt.plot(procs, speedup, marker='o', label=f"Pivot {pivot}")

        # 理想加速比
        if not ideal_line_drawn:
            plt.plot(procs, procs, linestyle='--', color='gray', label="Ideal Speedup")
            ideal_line_drawn = True

    plt.title(f"Speedup vs. Number of Processes\n{input_name}")
    plt.xlabel("Number of Processes")
    plt.ylabel("Speedup")
    plt.xticks(procs)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"speedup_{input_name.replace(' ', '_')}.png")
    plt.show()
