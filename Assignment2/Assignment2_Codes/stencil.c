#include "stencil.h"
#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define STENCIL_WIDTH 5
#define EXTENT (STENCIL_WIDTH/2)  // radius = 2

int read_input(const char *file_name, double **values);
int write_output(char *file_name, const double *output, int num_values);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s input_file output_file num_steps\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    // 1. 解析命令行参数
    char *input_file = argv[1];
    char *output_file = argv[2];
    int num_steps = atoi(argv[3]);

    // 2. 仅 rank 0 读取输入数据
    double *global_input = NULL;
    int num_values = 0;
    if (rank == 0) {
        if ((num_values = read_input(input_file, &global_input)) < 0) {
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
    }

    // 3. 广播 num_values
    MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 4. 设置与 Stencil 相关的参数
    double h = 2.0 * PI / num_values;
    // Centered difference stencil (width = 5)
    const double STENCIL[STENCIL_WIDTH] = {
        1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)
    };

    // 5. 计算本地数据大小
    int base = num_values / size;
    int rem  = num_values % size;
    int local_N = (rank < rem) ? (base + 1) : base;

    // 6. 计算 Scatter/Gather 的分发信息
    int *sendcounts = (int*)malloc(size * sizeof(int));
    int *displs     = (int*)malloc(size * sizeof(int));
    if (!sendcounts || !displs) {
        perror("malloc failed");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    int offset = 0;
    for (int i = 0; i < size; i++) {
        sendcounts[i] = (i < rem) ? (base + 1) : base;
        displs[i] = offset;
        offset += sendcounts[i];
    }

    // 7. 为 local_input 分配 (local_N + 2*EXTENT) 的空间，以容纳左右 ghost cells
    double *local_input = (double*)malloc((local_N + 2*EXTENT) * sizeof(double));
    if (!local_input) {
        perror("malloc failed");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    // 8. 为计算结果分配一个只含本地大小的临时缓冲
    //    （存放计算后结果，再复制到 local_input 中间）
    double *local_output = (double*)malloc(local_N * sizeof(double));
    if (!local_output) {
        perror("malloc failed");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    // 9. 将全局数据按对应分块分发至各个进程（中间区域是 [EXTENT, EXTENT+local_N)）
    MPI_Scatterv(
        global_input, sendcounts, displs, MPI_DOUBLE,
        &local_input[EXTENT],  // 从 EXTENT 开始存放
        local_N, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    // rank 0 不再需要 global_input
    if (rank == 0) {
        free(global_input);
        global_input = NULL;
    }

    // 10. 计算左右邻居的 rank（周期边界）
    int left_rank  = (rank - 1 + size) % size;
    int right_rank = (rank + 1) % size;

    // 11. 计时开始
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // 12. 主循环：进行 num_steps 次迭代
    for (int step = 0; step < num_steps; step++) {
        // 12.1 交换 ghost cells
        // 将 [EXTENT, EXTENT + EXTENT) 发给左邻居；接收右邻居的数据到 [local_N + EXTENT, local_N + 2*EXTENT)
        MPI_Sendrecv(&local_input[EXTENT],  // 发送起点
                     EXTENT, MPI_DOUBLE,    // 发送大小
                     left_rank, 0,          // 目标和 tag
                     &local_input[local_N + EXTENT],  // 接收起点
                     EXTENT, MPI_DOUBLE,    // 接收大小
                     right_rank, 0,         // 源和 tag
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // 将 [local_N, local_N + EXTENT) 发给右邻居；接收左邻居的数据到 [0, EXTENT)
        MPI_Sendrecv(&local_input[local_N], 
                     EXTENT, MPI_DOUBLE,
                     right_rank, 0,
                     &local_input[0], 
                     EXTENT, MPI_DOUBLE,
                     left_rank, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // 12.2 计算本地区域
        // local_input 的可用数据区间是 [0, local_N + 2*EXTENT)
        // 我们要在 [EXTENT, EXTENT + local_N) 上计算
        for (int i = 0; i < local_N; i++) {
            double val = 0.0;
            // i + j 访问时，要注意起点是 EXTENT
            // local_input[EXTENT + i + j - EXTENT], 
            // 实际就是 local_input[i + j].
            for (int j = 0; j < STENCIL_WIDTH; j++) {
                val += STENCIL[j] * local_input[i + j];
            }
            local_output[i] = val;
        }

        // 12.3 将计算结果复制回 local_input 的有效区域，以供下次迭代
        // 下次迭代的 ghost 交换仍基于 local_input，所以需要更新 [EXTENT, EXTENT + local_N)
        if (step < num_steps - 1) {
            memcpy(&local_input[EXTENT], local_output, local_N * sizeof(double));
        }
    }

    // 13. 最终结果收集（local_output 中还存着最后一次迭代的结果）
    // 但由于我们在最后一次迭代并未将 local_output 再次拷贝回 local_input，
    // 最后一轮计算结果就在 local_output 里。可任意选择其中一个做 Gatherv。
    // 这里我们就从 local_output 来收集更简单。
    double *final_output = NULL;
    if (rank == 0) {
        final_output = (double*)malloc(num_values * sizeof(double));
        if (!final_output) {
            perror("malloc failed");
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
    }

    MPI_Gatherv(
        local_output,        // 发送起点
        local_N, MPI_DOUBLE, // 发送大小
        final_output,        // 接收起点（rank 0 的全局数组）
        sendcounts, displs, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    // 14. 计时结束 & 输出结果
    double end_time = MPI_Wtime();
    double local_elapsed = end_time - start_time;
    double max_elapsed;
    MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Max Time across ranks: %.6f seconds\n", max_elapsed);

#ifdef PRODUCE_OUTPUT_FILE
        // 写文件
        if (write_output(output_file, final_output, num_values) != 0) {
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
#endif
        free(final_output);
    }

    // 15. 释放资源
    free(local_input);
    free(local_output);
    free(sendcounts);
    free(displs);

    MPI_Finalize();
    return 0;
}

/* ==========================
 * read_input / write_output
 * ========================== */

int read_input(const char *file_name, double **values) {
    FILE *file;
    if (NULL == (file = fopen(file_name, "r"))) {
        perror("Couldn't open input file");
        return -1;
    }
    int num_values;
    if (EOF == fscanf(file, "%d", &num_values)) {
        perror("Couldn't read element count from input file");
        fclose(file);
        return -1;
    }
    *values = (double*)malloc(num_values * sizeof(double));
    if (!(*values)) {
        perror("malloc for input data failed");
        fclose(file);
        return -1;
    }
    for (int i = 0; i < num_values; i++) {
        if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
            perror("Couldn't read elements from input file");
            fclose(file);
            return -1;
        }
    }
    if (0 != fclose(file)) {
        perror("Warning: couldn't close input file");
    }
    return num_values;
}

int write_output(char *file_name, const double *output, int num_values) {
    FILE *file;
    if (NULL == (file = fopen(file_name, "w"))) {
        perror("Couldn't open output file");
        return -1;
    }
    for (int i = 0; i < num_values; i++) {
        if (fprintf(file, "%.4f ", output[i]) < 0) {
            perror("Couldn't write data");
        }
    }
    if (fprintf(file, "\n") < 0) {
        perror("Couldn't write newline");
    }
    if (0 != fclose(file)) {
        perror("Warning: couldn't close output file");
    }
    return 0;
}
