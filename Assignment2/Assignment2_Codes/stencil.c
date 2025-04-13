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

    // 3. 广播 num_values 给所有进程
    MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 4. 设置 stencil 参数
    double h = 2.0 * PI / num_values;
    const double STENCIL[STENCIL_WIDTH] = {
        1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)
    };

    // 5. 计算每个进程的局部数据大小
    int base = num_values / size;
    int rem  = num_values % size;
    int local_N = (rank < rem) ? (base + 1) : base;

    // 6. 设置 Scatter/Gather 的分发信息（sendcounts 和 displs）
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

    // 7. 分配双缓冲区，每个大小为 (local_N + 2*EXTENT)
    double *buf_current = (double*)malloc((local_N + 2*EXTENT) * sizeof(double));
    double *buf_next    = (double*)malloc((local_N + 2*EXTENT) * sizeof(double));
    if (!buf_current || !buf_next) {
        perror("malloc failed");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    // 8. 将全局数据分发到各进程，
    //    有效数据存储在 buf_current[EXTENT, EXTENT+local_N)
    MPI_Scatterv(
        global_input, sendcounts, displs, MPI_DOUBLE,
        &buf_current[EXTENT],
        local_N, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );
    if (rank == 0) {
        free(global_input);
        global_input = NULL;
    }

    // 9. 计算左右邻居的 rank（周期边界）
    int left_rank  = (rank - 1 + size) % size;
    int right_rank = (rank + 1) % size;

    // 10. 定义缓冲区中有效数据的下标范围
    int eff_start = EXTENT;
    int eff_end   = EXTENT + local_N;  // 有效区长度 = local_N

    // 11. 确保初始数据分发完成后开始计时（计时不包含输入和初始分发）
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // 12. 主循环：进行 num_steps 次迭代计算
    for (int step = 0; step < num_steps; step++) {
        MPI_Request reqs[4];
        // 非阻塞接收更新 ghost 区域
        // 接收右侧 ghost：存放在 buf_current[eff_end, eff_end+EXTENT)
        MPI_Irecv(&buf_current[eff_end], EXTENT, MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD, &reqs[0]);
        // 接收左侧 ghost：存放在 buf_current[0, EXTENT)
        MPI_Irecv(&buf_current[0], EXTENT, MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD, &reqs[1]);

        // 非阻塞发送自身边界数据
        // 发送左边界：buf_current[eff_start, eff_start+EXTENT)
        MPI_Isend(&buf_current[eff_start], EXTENT, MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD, &reqs[2]);
        // 发送右边界：buf_current[eff_end-EXTENT, eff_end)
        MPI_Isend(&buf_current[eff_end - EXTENT], EXTENT, MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD, &reqs[3]);

        // 与通信重叠，先计算内部安全区域（不依赖最新 ghost 数据），
        // 即有效区 [eff_start+EXTENT, eff_end-EXTENT)
        for (int i = eff_start + EXTENT; i < eff_end - EXTENT; i++) {
            double sum = 0.0;
            for (int j = -EXTENT; j <= EXTENT; j++) {
                sum += STENCIL[j + EXTENT] * buf_current[i + j];
            }
            buf_next[i] = sum;
        }

        // 等待非阻塞通信完成后，再计算边界区域
        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

        // 计算左边界区域：有效区 [eff_start, eff_start+EXTENT)
        for (int i = eff_start; i < eff_start + EXTENT; i++) {
            double sum = 0.0;
            for (int j = -EXTENT; j <= EXTENT; j++) {
                sum += STENCIL[j + EXTENT] * buf_current[i + j];
            }
            buf_next[i] = sum;
        }
        // 计算右边界区域：有效区 [eff_end-EXTENT, eff_end)
        for (int i = eff_end - EXTENT; i < eff_end; i++) {
            double sum = 0.0;
            for (int j = -EXTENT; j <= EXTENT; j++) {
                sum += STENCIL[j + EXTENT] * buf_current[i + j];
            }
            buf_next[i] = sum;
        }

        // 利用双缓冲策略交换 buf_current 与 buf_next 指针，无需 memcpy
        double *tmp = buf_current;
        buf_current = buf_next;
        buf_next = tmp;
    }

    // 13. 迭代结束后，调用 Barrier 确保所有进程完成计算，
    //    然后计时结束（不包含最终数据收集时间）
    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    double local_elapsed = end_time - start_time;
    double max_elapsed;
    MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // 14. 在计时结束后收集最终结果（不计入计时）
    double *final_output = NULL;
    if (rank == 0) {
        final_output = (double*)malloc(num_values * sizeof(double));
        if (!final_output) {
            perror("malloc failed");
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
    }
    MPI_Gatherv(
        &buf_current[eff_start],  // 仅发送有效区域数据
        local_N, MPI_DOUBLE,
        final_output, sendcounts, displs, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    // 15. rank 0 输出计时结果，并在定义 PRODUCE_OUTPUT_FILE 时写入输出文件
    if (rank == 0) {
        printf("%f\n", max_elapsed);
#ifdef PRODUCE_OUTPUT_FILE
        if (write_output(output_file, final_output, num_values) != 0) {
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
#endif
        free(final_output);
    }

    // 16. 释放分配的资源
    free(buf_current);
    free(buf_next);
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
    if ((file = fopen(file_name, "r")) == NULL) {
        perror("Couldn't open input file");
        return -1;
    }
    int num_values;
    if (fscanf(file, "%d", &num_values) == EOF) {
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
        if (fscanf(file, "%lf", &((*values)[i])) == EOF) {
            perror("Couldn't read elements from input file");
            fclose(file);
            return -1;
        }
    }
    if (fclose(file) != 0) {
        perror("Warning: couldn't close input file");
    }
    return num_values;
}

int write_output(char *file_name, const double *output, int num_values) {
    FILE *file;
    if ((file = fopen(file_name, "w")) == NULL) {
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
    if (fclose(file) != 0) {
        perror("Warning: couldn't close output file");
    }
    return 0;
}
