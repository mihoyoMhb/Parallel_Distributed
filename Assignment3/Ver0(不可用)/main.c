// main.c
#include "quicksort.h"
#include "pivot.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // 参数检查
    if (argc != 4) {
        if (rank == ROOT)
            fprintf(stderr, "Usage: %s <input> <output> <pivot_strategy>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    
    char *input_file = argv[1];
    char *output_file = argv[2];
    int pivot_strategy = atoi(argv[3]);
    
    // 读取输入数据（仅在根进程）
    int *all_elements = NULL, *my_elements = NULL;
    int n = 0, local_n = 0;
    if (rank == ROOT) {
        n = read_input(input_file, &all_elements);
        if (n < 0) {
            fprintf(stderr, "Error reading input file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    
    // 分发数据
    local_n = distribute_from_root(all_elements, n, &my_elements);
    
    // 本地排序（步骤2）
    qsort(my_elements, local_n, sizeof(int), compare);
    
    // 测量排序时间
    double start_time = MPI_Wtime();
    
    // 全局排序（步骤3）
    local_n = global_sort(&my_elements, local_n, MPI_COMM_WORLD, pivot_strategy);
    
    double end_time = MPI_Wtime();
    
    // 收集结果到根进程
    if (rank == ROOT) {
        all_elements = realloc(all_elements, n * sizeof(int));
    }
    gather_on_root(all_elements, my_elements, local_n);
    
    // 验证和输出
    if (rank == ROOT) {
        check_and_print(all_elements, n, output_file);
        printf("%.6f\n", end_time - start_time);
    }
    
    // 清理资源
    free(my_elements);
    if (rank == ROOT) free(all_elements);
    
    MPI_Finalize();
    return 0;
}