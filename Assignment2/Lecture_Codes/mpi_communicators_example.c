// File: mpi_communicators_example.c
#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Comm new_comm;  // 新的communicator，用于存储奇数或偶数组

    // 初始化MPI环境
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // 获取在MPI_COMM_WORLD中的排名
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // 获取总进程数

    // 假设我们要求至少两个进程进行分组
    if (size < 2) {
        if (rank == 0)
            printf("请使用至少两个进程运行此示例。\n");
        MPI_Finalize();
        return 1;
    }

    // 根据进程的奇偶性设置color（奇数为1，偶数为0，或者可以写模N来分组）
    
    // int color = rank % N;  // N为分组数
    int color = rank % 2;  
    // key用于控制新communicator中进程的排序，此处直接使用rank
    int key = rank;
    // 利用MPI_Comm_split根据color分组，创建新communicator new_comm
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &new_comm);

    // 获取在新communicator中的新的排名和进程数
    int new_rank, new_size;
    MPI_Comm_rank(new_comm, &new_rank);
    MPI_Comm_size(new_comm, &new_size);

    // 每个进程赋予一个值（例如原rank+1），在各子组内执行求和归约
    int local_value = rank + 1;
    // 如果是用new_rank，那么结果不一样，new_rank是新comm中的rank: 0,1,2...
    //int local_value = new_rank + 1;
    int sum;
    MPI_Reduce(&local_value, &sum, 1, MPI_INT, MPI_SUM, 0, new_comm);

    // 在各自的子communicator的根进程打印归约结果
    if (new_rank == 0) {
        printf("在新communicator中，color = %d 的组求和结果为：%d (进程组总数：%d)\n", color, sum, new_size);
    }

    // 释放新创建的communicator
    MPI_Comm_free(&new_comm);
    MPI_Finalize();
    return 0;
}
