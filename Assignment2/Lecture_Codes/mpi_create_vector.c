#include <mpi.h>
#include <stdio.h>

#define N 4  // 矩阵维度

int main(int argc, char *argv[]) {
    int rank, size;
    int matrix[N][N];
    MPI_Datatype diag_type;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 仅允许至少两个进程进行演示
    if (size < 2) {
        if (rank == 0)
            printf("请使用至少两个进程运行此示例.\n");
        MPI_Finalize();
        return 1;
    }

    // 假设rank 0填充数据
    if (rank == 0) {
        int counter = 1;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                matrix[i][j] = counter++;
            }
        }

        // 打印矩阵，便于验证
        printf("Matrix at Rank 0:\n");
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                printf("%3d ", matrix[i][j]);
            }
            printf("\n");
        }
    }

    // 创建派生数据类型，用于选择对角线元素
    // count = N，即有N个对角线元素
    // blocklength = 1，每个块只有一个元素
    // stride = N+1，因为从matrix[0][0]到matrix[1][1]需要跳过N - 1个非对角元素和跨过一个位置
    MPI_Type_vector(N, 1, N + 1, MPI_INT, &diag_type);
    MPI_Type_commit(&diag_type);

    if (rank == 0) {
        // 将矩阵的对角线元素发送给Rank 1
        MPI_Send(&matrix[0][0], 1, diag_type, 1, 0, MPI_COMM_WORLD);
        printf("Rank 0 sent diagonal elements to Rank 1\n");
    } else if (rank == 1) {
        int diag_elements[N];
        // 接收数据到一维数组，注意个数是N个整数
        MPI_Recv(diag_elements, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Rank 1 received diagonal elements: ");
        for (int i = 0; i < N; i++) {
            printf("%d ", diag_elements[i]);
        }
        printf("\n");
    }

    MPI_Type_free(&diag_type);
    MPI_Finalize();
    return 0;
}