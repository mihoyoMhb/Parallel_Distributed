#ifndef QUICKSORT_H_
#define QUICKSORT_H_

#include <mpi.h>

/**
 * Verify that elements are sorted in ascending order. If not, write an error
 * message to stdout. Thereafter, print all elements (in the order they are
 * stored) to a file named according to the last argument.
 * @param elements Elements to check and print
 * @param n Number of elements
 * @param file_name Name of output file
 * @return 0 on success, -2 on I/O error
 */
int check_and_print(int *elements, int n, char *file_name);

/**
 * Distribute all elements from root to the other processes as evenly as
 * possible. Note that this method allocates memory for my_elements. This must
 * be freed by the caller!
 * @param all_elements Elements to distribute (Not significant in other processes)
 * @param n Number of elements in all_elements
 * @param my_elements Pointer to buffer where the local elements will be stored
 * @return Number of elements received by the current process
 */
int distribute_from_root(int *all_elements, int n, int **my_elements);

/**
 * Gather elements from all processes on root. Put root's elements first and
 * thereafter elements from the other nodes in the order of their ranks (so that
 * elements from process i come after the elements from process i-1).
 * @param all_elements Buffer on root where the elements will be stored
 * @param my_elements Elements to be gathered from the current process
 * @param local_n Number of elements in my_elements
 */
void gather_on_root(int *all_elements, int *my_elements, int local_n);

/**
 * Perform the global part of parallel quick sort. This function assumes that
 * the elements is sorted within each node. When the function returns, all
 * elements owned by process i are smaller than or equal to all elements owned
 * by process i+1, and the elements are sorted within each node.
 * @param elements Pointer to the array of sorted values on the current node. Will point to a(n) (new) array with the sorted elements when the function returns.
 * @param n Length of *elements
 * @param MPI_Comm Communicator containing all processes participating in the global sort
 * @param pivot_strategy Tells how to select the pivot element. See documentation of select_pivot in pivot.h.
 * @return New length of *elements
 */
int global_sort(int **elements, int n, MPI_Comm, int pivot_strategy);

/**
 * Merge v1 and v2 to one array, sorted in ascending order, and store the result
 * in result.
 * @param v1 Array to merge
 * @param n1 Length of v1
 * @param v2 Array to merge
 * @param n2 Length of v2
 * @param result Array for merged result (must be allocated before!)
 */
void merge_ascending(int *v1, int n1, int *v2, int n2, int *result);

/**
 * Read problem size and elements from the file whose name is given as an
 * argument to the function and populate elements accordingly. Note that this
 * method allocates memory for elements. This must be freed by the caller!
 * @param file_name Name of input file
 * @param elements Pointer to array to be created and populated
 * @return Number of elements read and stored
 */
int read_input(char *file_name, int **elements);

/**
 * Check if a number of elements are sorted in ascending order. If they aren't,
 * print an error message specifying the first two elements that are in wrong
 * order.
 * @param elements Array to check
 * @param n Length of elements
 * @return 1 if elements is sorted in ascending order, 0 otherwise
 */
int sorted_ascending(int *elements, int n);

/**
 * Swap the values pointed at by e1 and e2.
 */
void swap(int *e1, int *e2);


// /**
//  * 选择轴心元素并确定局部小于等于轴心和大于轴心的元素数量
//  * 
//  * @param elements_ptr 指向元素数组的指针
//  * @param n 元素数量
//  * @param communicator MPI通信器
//  * @param pivot_strategy 轴心选择策略
//  * @param num_small_local 返回小于等于轴心的元素数量
//  * @param num_large_local 返回大于轴心的元素数量
//  * @return 返回小于等于轴心的元素数量(k)
//  */
// int select_pivot_and_partition(int **elements_ptr, int n, MPI_Comm communicator, 
//                             int pivot_strategy, int *num_small_local, int *num_large_local);

                    
// /**
//  * 确定进程分组和通信伙伴
//  * 
//  * @param my_rank 当前进程的序号
//  * @param num_procs 总进程数
//  * @param num_small_local 小于等于轴心的元素数量
//  * @param num_large_local 大于轴心的元素数量
//  * @param elements_ptr 指向元素数组的指针
//  * @param color 返回进程的颜色(分组)
//  * @param partner_rank 返回通信伙伴的进程序号
//  * @param send_buf_ptr 返回要发送的缓冲区的起始位置
//  * @param send_count 返回要发送的元素数量
//  * @param keep_buf_start_ptr 返回要保留的缓冲区的起始位置
//  * @param keep_count 返回要保留的元素数量
//  */
// void determine_process_groups(int my_rank, int num_procs, int num_small_local, int num_large_local,
//                              int **elements_ptr, int *color, int *partner_rank,
//                              int **send_buf_ptr, int *send_count,
//                              int **keep_buf_start_ptr, int *keep_count);



// /**
//  * 交换数据：向通信伙伴发送数据并从伙伴接收数据
//  * 
//  * @param send_buf_ptr 发送缓冲区的起始位置
//  * @param send_count 发送的元素数量
//  * @param partner_rank 通信伙伴的进程序号
//  * @param communicator MPI通信器
//  * @param my_rank 当前进程的序号
//  * @param recv_count_ptr 返回接收到的元素数量
//  * @return 返回接收到的数据缓冲区，如无数据则为NULL
//  */
// int* exchange_data(int *send_buf_ptr, int send_count, int partner_rank, 
//                  MPI_Comm communicator, int my_rank, int *recv_count_ptr);

// /**
//  * 合并保留的数据和接收到的数据
//  * 
//  * @param keep_buf_start_ptr 保留缓冲区的起始位置
//  * @param keep_count 保留的元素数量
//  * @param recv_buf 接收到的缓冲区
//  * @param recv_count 接收到的元素数量
//  * @param my_rank 当前进程的序号
//  * @param communicator MPI通信器
//  * @param new_n_ptr 返回合并后的元素总数
//  * @return 返回合并后的新元素数组
//  */
// int* merge_data(int *keep_buf_start_ptr, int keep_count, int *recv_buf, int recv_count,
//                int my_rank, MPI_Comm communicator, int *new_n_ptr);
               

// /**
//  * 创建新的通信器并进行递归排序
//  * 
//  * @param elements_ptr 指向元素数组的指针
//  * @param n 元素数量
//  * @param communicator 当前MPI通信器
//  * @param color 进程颜色(分组)
//  * @param my_rank 当前进程的序号
//  * @param pivot_strategy 轴心选择策略
//  * @return 返回排序后的元素数量
//  */
// int recursive_sort(int **elements_ptr, int n, MPI_Comm communicator, 
//                   int color, int my_rank, int pivot_strategy);
#endif