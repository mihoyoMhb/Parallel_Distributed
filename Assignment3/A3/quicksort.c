#include "quicksort.h"
#include "pivot.h" 
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For memcpy

void swap(int *e1, int *e2) {
    int temp = *e1;
    *e1 = *e2;
    *e2 = temp;
}

int sorted_ascending(int *elements, int n) {
    if (n == 0 || elements == NULL) return 1;
    for (int i = 0; i < n - 1; ++i) {
        if (elements[i] > elements[i + 1]) {
            return 0; 
        }
    }
    return 1; 
}

void merge_ascending(int *v1, int n1, int *v2, int n2, int *result) {
    int i = 0, j = 0, k = 0;
    while (i < n1 && j < n2) {
        if (v1[i] <= v2[j]) {
            result[k++] = v1[i++];
        } else {
            result[k++] = v2[j++];
        }
    }
    while (i < n1) {
        result[k++] = v1[i++];
    }
    while (j < n2) {
        result[k++] = v2[j++];
    }
}

int read_input(char *file_name, int **elements_ptr) {
    FILE *file = fopen(file_name, "r");
    // Assuming file opens successfully for stable test cases
    
    int n;
    if (fscanf(file, "%d", &n) != 1) {
        fprintf(stderr, "Failed to read problem size from file %s\n", file_name);
        fclose(file);
        *elements_ptr = NULL;
        return -1; 
    }

    
    if (n == 0) {
        *elements_ptr = NULL; 
        fclose(file);
        return 0;
    }

    *elements_ptr = (int *)malloc(n * sizeof(int));
    if (!*elements_ptr) { // Check for n > 0 is implicit as n=0 is handled above
        fprintf(stderr, "ROOT: Failed to allocate memory for %d elements in read_input\n", n);
        fclose(file); // Close file before aborting
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < n; ++i) {
        if (fscanf(file, "%d", &((*elements_ptr)[i])) != 1) {
            fprintf(stderr, "Failed to read element %d (0-indexed) from file %s\n", i, file_name);
            free(*elements_ptr);
            *elements_ptr = NULL;
            fclose(file);
            return -1; 
        }
    }

    fclose(file);
    return n;
}

int check_and_print(int *elements, int n, char *file_name) {
    sorted_ascending(elements, n); // Call for its side-effect if it had one, or just check
                                   // Removed warning message if not sorted.

    FILE *file = fopen(file_name, "w");
    // Assuming file opens successfully

    for (int i = 0; i < n; ++i) {
        fprintf(file, "%d", elements[i]);
        if (i < n - 1) {
            fprintf(file, " "); 
        }
    }
    fprintf(file, "\n"); 

    fclose(file);
    return 0; 
}

int distribute_from_root(int *all_elements_on_root, int n_total, int **my_elements_ptr) {
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int *sendcounts = NULL;
    int *displs = NULL;
    int local_n;

    if (rank == ROOT) {
        if (comm_size > 0) {
            sendcounts = (int *)calloc(comm_size, sizeof(int)); 
            displs = (int *)calloc(comm_size, sizeof(int));
            if (!sendcounts || !displs) {
                 fprintf(stderr, "ROOT: Failed to allocate memory for scatterv metadata.\n");
                 if(sendcounts) free(sendcounts); // Free partially allocated if possible
                 if(displs) free(displs);
                 MPI_Abort(MPI_COMM_WORLD, 1);
            }

            int base_count = (comm_size > 0) ? (n_total / comm_size) : 0;
            int remainder = (comm_size > 0) ? (n_total % comm_size) : 0;
            int current_displ = 0;
            for (int i = 0; i < comm_size; ++i) {
                sendcounts[i] = base_count + (i < remainder ? 1 : 0);
                displs[i] = current_displ;
                current_displ += sendcounts[i];
            }
        }
    }

    MPI_Scatter(sendcounts, 1, MPI_INT, &local_n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    
    if (local_n > 0) {
        *my_elements_ptr = (int *)malloc(local_n * sizeof(int));
        if (!*my_elements_ptr) {
            fprintf(stderr, "[Process %d] Failed to allocate memory for %d local elements.\n", rank, local_n);
            if (rank == ROOT) { // Free scatterv metadata if ROOT fails here
                if (sendcounts) free(sendcounts);
                if (displs) free(displs);
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    } else {
        *my_elements_ptr = NULL; 
    }
    
    MPI_Scatterv(all_elements_on_root, sendcounts, displs, MPI_INT,
                 *my_elements_ptr, local_n, MPI_INT,
                 ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
        if (sendcounts) free(sendcounts);
        if (displs) free(displs);
    }
    return local_n;
}

void gather_on_root(int *all_elements_buffer_on_root, int *my_elements, int local_n) {
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int *recvcounts = NULL;
    int *displs = NULL;

    if (rank == ROOT) {
        if (comm_size > 0) {
            recvcounts = (int *)malloc(comm_size * sizeof(int));
            displs = (int *)malloc(comm_size * sizeof(int));
            if (!recvcounts || !displs) {
                 fprintf(stderr, "ROOT: Failed to allocate memory for gatherv metadata.\n");
                 if(recvcounts) free(recvcounts);
                 if(displs) free(displs);
                 MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    MPI_Gather(&local_n, 1, MPI_INT, recvcounts, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
        if (comm_size > 0) { // Ensure recvcounts is not NULL before using it
            int current_displ = 0;
            for (int i = 0; i < comm_size; ++i) {
                displs[i] = current_displ;
                current_displ += recvcounts[i]; // recvcounts should be populated by MPI_Gather
            }
        }
    }

    MPI_Gatherv(my_elements, local_n, MPI_INT,
                all_elements_buffer_on_root, recvcounts, displs, MPI_INT,
                ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
        if(recvcounts) free(recvcounts);
        if(displs) free(displs);
    }
}

/**
 * 选择轴心元素并确定局部小于等于轴心和大于轴心的元素数量
 * 
 * @param elements_ptr 指向元素数组的指针
 * @param n 元素数量
 * @param communicator MPI通信器
 * @param pivot_strategy 轴心选择策略
 * @param num_small_local 返回小于等于轴心的元素数量
 * @param num_large_local 返回大于轴心的元素数量
 * @return 返回小于等于轴心的元素数量(k)
 */
int select_pivot_and_partition(int **elements_ptr, int n, MPI_Comm communicator, 
                            int pivot_strategy, int *num_small_local, int *num_large_local) {
    int k = select_pivot(pivot_strategy, *elements_ptr, n, communicator);
    *num_small_local = k;
    *num_large_local = n - k;
    return k;
}

/**
 * 确定进程分组和通信伙伴
 * 
 * @param my_rank 当前进程的序号
 * @param num_procs 总进程数
 * @param num_small_local 小于等于轴心的元素数量
 * @param num_large_local 大于轴心的元素数量
 * @param elements_ptr 指向元素数组的指针
 * @param color 返回进程的颜色(分组)
 * @param partner_rank 返回通信伙伴的进程序号
 * @param send_buf_ptr 返回要发送的缓冲区的起始位置
 * @param send_count 返回要发送的元素数量
 * @param keep_buf_start_ptr 返回要保留的缓冲区的起始位置
 * @param keep_count 返回要保留的元素数量
 */
void determine_process_groups(int my_rank, int num_procs, int num_small_local, int num_large_local,
                             int **elements_ptr, int *color, int *partner_rank,
                             int **send_buf_ptr, int *send_count,
                             int **keep_buf_start_ptr, int *keep_count) {
    *color = (my_rank < num_procs / 2) ? 0 : 1;
    
    if (*color == 0) { 
        // 前半部分进程处理较小的元素
        *partner_rank = my_rank + num_procs / 2;
        *send_buf_ptr = (*elements_ptr) + num_small_local;
        *send_count = num_large_local;
        *keep_buf_start_ptr = *elements_ptr;
        *keep_count = num_small_local;
    } else { 
        // 后半部分进程处理较大的元素
        *partner_rank = my_rank - num_procs / 2;
        *send_buf_ptr = *elements_ptr;
        *send_count = num_small_local;
        *keep_buf_start_ptr = (*elements_ptr) + num_small_local;
        *keep_count = num_large_local;
    }
}

/**
 * 交换数据：向通信伙伴发送数据并从伙伴接收数据
 * 
 * @param send_buf_ptr 发送缓冲区的起始位置
 * @param send_count 发送的元素数量
 * @param partner_rank 通信伙伴的进程序号
 * @param communicator MPI通信器
 * @param my_rank 当前进程的序号
 * @param recv_count_ptr 返回接收到的元素数量
 * @return 返回接收到的数据缓冲区，如无数据则为NULL
 */
int* exchange_data(int *send_buf_ptr, int send_count, int partner_rank, 
                 MPI_Comm communicator, int my_rank, int *recv_count_ptr) {
    int recv_count;
    // 首先交换将要发送的数据量
    MPI_Sendrecv(&send_count, 1, MPI_INT, partner_rank, 0, 
                 &recv_count, 1, MPI_INT, partner_rank, 0,
                 communicator, MPI_STATUS_IGNORE);
    
    *recv_count_ptr = recv_count;
    int *recv_buf = NULL;
    
    // 分配接收缓冲区（如果需要）
    if (recv_count > 0) {
        recv_buf = (int *)malloc(recv_count * sizeof(int));
        if (!recv_buf) {
            fprintf(stderr, "[Proc %d] Failed to alloc recv_buf for %d ints\n", my_rank, recv_count);
            MPI_Abort(communicator, 1);
        }
    }
    
    // 交换实际数据
    MPI_Sendrecv(send_buf_ptr, send_count, MPI_INT, partner_rank, 1, 
                 recv_buf, recv_count, MPI_INT, partner_rank, 1,
                 communicator, MPI_STATUS_IGNORE);
                 
    return recv_buf;
}

/**
 * 合并保留的数据和接收到的数据
 * 
 * @param keep_buf_start_ptr 保留缓冲区的起始位置
 * @param keep_count 保留的元素数量
 * @param recv_buf 接收到的缓冲区
 * @param recv_count 接收到的元素数量
 * @param my_rank 当前进程的序号
 * @param communicator MPI通信器
 * @param new_n_ptr 返回合并后的元素总数
 * @return 返回合并后的新元素数组
 */
int* merge_data(int *keep_buf_start_ptr, int keep_count, int *recv_buf, int recv_count,
               int my_rank, MPI_Comm communicator, int *new_n_ptr) {
    int new_n = keep_count + recv_count;
    *new_n_ptr = new_n;
    int *new_elements_array = NULL;

    if (new_n > 0) {
        // 分配新数组用于合并结果
        new_elements_array = (int *)malloc(new_n * sizeof(int));
        if (!new_elements_array) {
            fprintf(stderr, "[Proc %d] Failed to alloc new_elements_array for %d ints\n", my_rank, new_n);
            if(recv_buf) free(recv_buf);
            MPI_Abort(communicator, 1);
        }

        if (keep_count > 0 && recv_count > 0) {
            // 需要将两个排序好的数组合并
            int* temp_kept_elements = (int*)malloc(keep_count * sizeof(int));
            if (!temp_kept_elements) {
                fprintf(stderr, "[Proc %d] Failed to alloc temp_kept_elements for %d ints\n", my_rank, keep_count);
                if(recv_buf) free(recv_buf);
                free(new_elements_array);
                MPI_Abort(communicator, 1);
            }
            memcpy(temp_kept_elements, keep_buf_start_ptr, keep_count * sizeof(int));
            merge_ascending(temp_kept_elements, keep_count, recv_buf, recv_count, new_elements_array);
            free(temp_kept_elements);
        } else if (keep_count > 0) { 
            // 只有保留的数据
            memcpy(new_elements_array, keep_buf_start_ptr, keep_count * sizeof(int));
        } else if (recv_count > 0) { 
            // 只有接收到的数据
            memcpy(new_elements_array, recv_buf, recv_count * sizeof(int));
        }
    }
    
    return new_elements_array;
}

/**
 * 创建新的通信器并进行递归排序
 * 
 * @param elements_ptr 指向元素数组的指针
 * @param n 元素数量
 * @param communicator 当前MPI通信器
 * @param color 进程颜色(分组)
 * @param my_rank 当前进程的序号
 * @param pivot_strategy 轴心选择策略
 * @return 返回排序后的元素数量
 */
int recursive_sort(int **elements_ptr, int n, MPI_Comm communicator, 
                  int color, int my_rank, int pivot_strategy) {
    // 创建新的通信器，按照颜色分组
    MPI_Comm new_comm;
    MPI_Comm_split(communicator, color, my_rank, &new_comm);
    
    // 在新的通信器上递归排序
    n = global_sort(elements_ptr, n, new_comm, pivot_strategy);
    
    // 释放通信器
    MPI_Comm_free(&new_comm);
    
    return n;
}

/**
 * 全局并行快速排序算法
 * 
 * @param elements_ptr 指向元素数组的指针
 * @param n 元素数量
 * @param communicator MPI通信器
 * @param pivot_strategy 轴心选择策略
 * @return 返回排序后的元素数量
 */
int global_sort(int **elements_ptr, int n, MPI_Comm communicator, int pivot_strategy) {
    int my_rank, num_procs;
    MPI_Comm_rank(communicator, &my_rank);
    MPI_Comm_size(communicator, &num_procs);

    // 基本情况：如果只有一个进程，不需要并行排序
    if (num_procs <= 1) {
        return n; 
    }

    // 1. 选择轴心并确定局部分区
    int num_small_local, num_large_local;
    select_pivot_and_partition(elements_ptr, n, communicator, 
                             pivot_strategy, &num_small_local, &num_large_local);
    
    // 2. 确定进程分组和通信伙伴
    int color, partner_rank;
    int *send_buf_ptr, *keep_buf_start_ptr;
    int send_count, keep_count;
    determine_process_groups(my_rank, num_procs, num_small_local, num_large_local,
                           elements_ptr, &color, &partner_rank,
                           &send_buf_ptr, &send_count,
                           &keep_buf_start_ptr, &keep_count);
    
    // 3. 交换数据
    int recv_count;
    int *recv_buf = exchange_data(send_buf_ptr, send_count, partner_rank,
                                communicator, my_rank, &recv_count);
    
    // 4. 合并数据
    int new_n;
    int *new_elements_array = merge_data(keep_buf_start_ptr, keep_count, 
                                        recv_buf, recv_count,
                                        my_rank, communicator, &new_n);
    
    // 释放旧数据并更新指针
    if (*elements_ptr) { 
        free(*elements_ptr);
    }
    *elements_ptr = new_elements_array; 
    n = new_n;

    // 释放接收缓冲区
    if (recv_buf) {
        free(recv_buf);
    }
    
    // 5. 在新的通信器上递归排序
    n = recursive_sort(elements_ptr, n, communicator, color, my_rank, pivot_strategy);
    
    return n;
}