
// File: mpi_create_struct.c

#include <mpi.h>
#include <stdio.h>
#include <stddef.h>  // 包含offsetof宏需要的头文件

// 定义一个结构体类型
typedef struct {
    int id;
    double value;
} MyData;

int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 构造MPI派生数据类型与 MyData 对应

    // 定义结构体的成员信息:
    // blocklengths：每个成员包含的元素数目（这里每个成员只有一个元素）
    int blocklengths[2] = {1, 1};
    // offsets：各成员在结构体中的偏移量，使用标准宏offsetof计算偏移
    MPI_Aint offsets[2];
    offsets[0] = offsetof(MyData, id);
    offsets[1] = offsetof(MyData, value);
    // types：各成员对应的MPI数据类型, 长度为2，分别对应int和double
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};

    // 创建MPI派生数据类型，与MyData结构体内存布局一致
    MPI_Datatype mpi_mydata_type;
    MPI_Type_create_struct(2, blocklengths, offsets, types, &mpi_mydata_type);
    // 提交数据类型，使其可用于通信
    MPI_Type_commit(&mpi_mydata_type);

    // 定义两个MyData类型的变量用于发送和接收
    MyData data_send, data_recv;
    data_send.id = rank;         // 例如，进程号设置为id
    data_send.value = rank * 1.0;  // 设置value为进程号的浮点值

    // 这里只让进程0和进程1进行通信
    if (size < 2) {
        if (rank == 0)
            printf("请使用至少两个进程运行此示例.\n");
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        // 进程0接收来自进程1的数据
        MPI_Recv(&data_recv, 1, mpi_mydata_type, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process 0 received MyData: id = %d, value = %.2f\n", data_recv.id, data_recv.value);
    } else if (rank == 1) {
        // 进程1发送数据给进程0
        MPI_Send(&data_send, 1, mpi_mydata_type, 0, 0, MPI_COMM_WORLD);
        printf("Process 1 sent MyData: id = %d, value = %.2f\n", data_send.id, data_send.value);
    }

    // 释放MPI创建的数据类型
    MPI_Type_free(&mpi_mydata_type);
    MPI_Finalize();
    return 0;
}
