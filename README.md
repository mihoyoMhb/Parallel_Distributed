# Parallel_Distributed
# MPI 函数速查完整文档

## 1. 初始化与终止

### `MPI_Init`
```c
int MPI_Init(int *argc, char ***argv)
```
功能: 初始化 MPI 环境

参数:

- `argc`: 命令行参数数量的指针
- `argv`: 命令行参数数组的指针

### `MPI_Finalize`
```c
int MPI_Finalize()
```
功能: 终止 MPI 环境

---

## 2. 进程信息

### `MPI_Comm_size`
```c
int MPI_Comm_size(MPI_Comm comm, int *size)
```
功能: 获取通信器中的进程总数

参数:

- `comm`: 通信器（如 MPI_COMM_WORLD）
- `size`: 输出进程总数

### `MPI_Comm_rank`
```c
int MPI_Comm_rank(MPI_Comm comm, int *rank)
```
功能: 获取当前进程的 ID

参数:

- `comm`: 通信器
- `rank`: 输出进程 ID (0 到 size-1)

---

## 3. 点对点通信

### `MPI_Send`
```c
int MPI_Send(void *buf, int count, MPI_Datatype datatype, 
             int dest, int tag, MPI_Comm comm)
```
功能: 阻塞发送数据

参数:

- `buf`: 发送数据缓冲区
- `count`: 元素数量
- `datatype`: 数据类型（如 MPI_INT）
- `dest`: 目标进程 ID
- `tag`: 消息标签
- `comm`: 通信器

### `MPI_Recv`
```c
int MPI_Recv(void *buf, int count, MPI_Datatype datatype,
             int source, int tag, MPI_Comm comm, MPI_Status *status)
```
功能: 阻塞接收数据

参数:

- `buf`: 接收数据缓冲区
- `count`: 最大元素数量
- `datatype`: 数据类型
- `source`: 源进程 ID（可使用 MPI_ANY_SOURCE）
- `tag`: 消息标签（可使用 MPI_ANY_TAG）
- `comm`: 通信器
- `status`: 状态对象（存储实际接收信息）

---

## 4. 集合通信

### 基本 Scatter/Gather

#### `MPI_Scatter`
```c
int MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, int recvcount, MPI_Datatype recvtype,
                int root, MPI_Comm comm)
```
功能: 从根进程分发等量数据块

参数:

- `sendbuf`: 根进程的发送缓冲区（其他进程忽略）
- `sendcount`: 每个进程接收的元素数
- `sendtype`: 发送数据类型
- `recvbuf`: 接收缓冲区
- `recvcount`: 接收元素数
- `recvtype`: 接收数据类型
- `root`: 根进程 ID
- `comm`: 通信器

#### `MPI_Gather`
```c
int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype,
               int root, MPI_Comm comm)
```
功能: 将数据收集到根进程

参数: 同 MPI_Scatter，方向相反

---

### 变体 Scatterv/Gatherv

#### `MPI_Scatterv`
```c
int MPI_Scatterv(void *sendbuf, int *sendcounts, int *displs,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount,
                 MPI_Datatype recvtype, int root, MPI_Comm comm)
```
功能: 分发不等量数据块

参数:

- `sendcounts`: 数组，指定每个进程接收的元素数
- `displs`: 数组，指定根进程中各进程数据的偏移量
- 其他参数同 MPI_Scatter

#### `MPI_Gatherv`
```c
int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, int *recvcounts, int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm)
```
功能: 收集不等量数据块

参数: 同 MPI_Scatterv，方向相反

---

### 广播与规约

#### `MPI_Bcast`
```c
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
              int root, MPI_Comm comm)
```
功能: 从根进程广播数据到所有进程

参数:

- `buffer`: 根进程的发送缓冲区/其他进程的接收缓冲区
- 其他参数同 MPI_Send

#### `MPI_Reduce`
```c
int MPI_Reduce(void *sendbuf, void *recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, int root,
               MPI_Comm comm)
```
功能: 执行规约操作（如求和、最大值）

参数:

- `op`: 规约操作（如 MPI_SUM, MPI_MAX）
- 其他参数同 MPI_Gather

---

## 5. 同步

### `MPI_Barrier`
```c
int MPI_Barrier(MPI_Comm comm)
```
功能: 阻塞直到所有进程到达此调用

---

## 6. 数据类型

常用预定义类型：

- `MPI_INT`
- `MPI_DOUBLE`
- `MPI_CHAR`
- `MPI_BYTE`

---

## 7. 关键注意事项

- **通信器一致性**：所有进程必须使用相同的通信器  
- **标签匹配**：发送/接收操作的标签必须一致  
- **死锁风险**：避免阻塞调用的循环依赖  
- **集合操作参与**：所有进程必须参与集合通信调用  
- **缓冲区安全**：通信过程中勿修改缓冲区  

**Scatterv/Gatherv 专用规则：**

- `sendcounts` 和 `displs` 只需根进程有效  
- 位移以元素为单位（非字节）  
- 接收方 `recvcount` 必须等于根进程 `sendcounts[rank]`  

---

## 8. 示例代码

### 基础模式
```c
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // MPI 操作...
    
    MPI_Finalize();
    return 0;
}
```

### Scatterv 示例
```c
int main() {
    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int *sendbuf = NULL;
    int *sendcounts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));
    
    if (rank == 0) { // 根进程准备数据
        int total = 0;
        for (int i=0; i<size; i++) {
            sendcounts[i] = i+1;  // 进程 i 接收 i+1 个元素
            total += sendcounts[i];
        }
        sendbuf = malloc(total * sizeof(int));
        for (int i=0; i<total; i++) sendbuf[i] = i;
        
        displs[0] = 0;
        for (int i=1; i<size; i++) {
            displs[i] = displs[i-1] + sendcounts[i-1];
        }
    }

    int recvcount;
    MPI_Scatter(sendcounts, 1, MPI_INT, &recvcount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *recvbuf = malloc(recvcount * sizeof(int));
    MPI_Scatterv(sendbuf, sendcounts, displs, MPI_INT,
                 recvbuf, recvcount, MPI_INT, 0, MPI_COMM_WORLD);

    free(sendcounts); free(displs); free(recvbuf);
    if (rank == 0) free(sendbuf);
    MPI_Finalize();
}
```

---

## 9. 使用场景

**Scatterv/Gatherv**:

- 不规则数据分布（如矩阵不等长行）
- 动态负载均衡
- 根进程预处理不同大小的数据块

**练习建议：** 尝试实现主从模式、域分解算法，并测试不同进程数下的行为。
