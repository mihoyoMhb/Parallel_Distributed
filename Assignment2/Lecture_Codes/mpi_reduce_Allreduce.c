#include <mpi.h>
#include <stdio.h>

// 定义一个结构体类型
typedef struct {
    int id;
    double value;
} MyData;

int main(int argc, char *argv[]) {
    int rank, size;
    // 待归约的变量
    double local_val, sum_val;
    double all_sum_val;

    // 初始化MPI环境
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 每个进程初始化不同的值
    local_val = rank + 1; // 例如：进程0为1，进程1为2，依此类推

    // 使用MPI_Reduce将所有local_val的值求和，并将结果发送到进程0
    MPI_Reduce(&local_val, &sum_val, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("MPI_Reduce: The sum of local values is %f\n", sum_val);
    }

    // 使用MPI_Allreduce将所有local_val的值求和，并让所有进程获得结果
    MPI_Allreduce(&local_val, &all_sum_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    printf("Rank %d: MPI_Allreduce result is %.0f\n", rank, all_sum_val);

    // 示例：利用MPI_Type_create_struct定义一个MPI派生数据类型用于结构体MyData

    // 定义结构体的成员信息
    int blocklengths[2] = {1, 1};
    MPI_Aint offsets[2];
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};

    // 计算成员相对于结构体开头的偏移量
    offsets[0] = offsetof(MyData, id);
    offsets[1] = offsetof(MyData, value);

    MPI_Datatype mpi_mydata_type;
    MPI_Type_create_struct(2, blocklengths, offsets, types, &mpi_mydata_type);
    MPI_Type_commit(&mpi_mydata_type);

    // 示例：构造一个MyData类型的变量，进行数据传输
    MyData data_send, data_recv;
    data_send.id = rank;
    data_send.value = rank * 1.0;  // 举例赋值

    if (rank == 0) {
        // 进程0接收来自进程1的数据
        MPI_Recv(&data_recv, 1, mpi_mydata_type, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process 0 received MyData: id = %d, value = %.2f\n", data_recv.id, data_recv.value);
    } else if (rank == 1) {
        // 进程1发送数据给进程0
        MPI_Send(&data_send, 1, mpi_mydata_type, 0, 0, MPI_COMM_WORLD);
    }

    // 释放MPI类型
    MPI_Type_free(&mpi_mydata_type);
    MPI_Finalize();
    return 0;
}