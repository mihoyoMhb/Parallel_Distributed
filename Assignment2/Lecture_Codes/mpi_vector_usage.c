#include <mpi.h>
#include <stdio.h>

#define ROWS 4
#define COLS 5

int main(int argc, char* argv[]) {
    int rank, size;
    int matrix[ROWS][COLS];
    MPI_Datatype column_type;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 假设rank 0初始化数据
    if (rank == 0) {
        for (int i = 0; i < ROWS; i++) {
            for (int j = 0; j < COLS; j++) {
                matrix[i][j] = i * COLS + j;
            }
        }
    }

    // 创建一个派生数据类型，表示二维数组中某一列的数据
    // 参数解释：count = ROWS (块数), blocklength = 1 (每块一个元素), 
    // stride = COLS (每行的元素个数，即跳过整行数据)，数据类型为MPI_INT
    MPI_Type_vector(ROWS, 1, COLS, MPI_INT, &column_type);
    MPI_Type_commit(&column_type);

    if (rank == 0) {
        // 发送二维数组中第2列的数据（例如索引为1的列）给rank 1
        MPI_Send(&matrix[0][1], 1, column_type, 1, 0, MPI_COMM_WORLD);
        printf("Rank 0 sent column data to Rank 1\n");
    } else if (rank == 1) {
        int recv_buffer[ROWS];
        MPI_Recv(recv_buffer, ROWS, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Rank 1 received column data: ");
        for (int i = 0; i < ROWS; i++) {
            printf("%d ", recv_buffer[i]);
        }
        printf("\n");
    }

    MPI_Type_free(&column_type);
    MPI_Finalize();

    return 0;
}