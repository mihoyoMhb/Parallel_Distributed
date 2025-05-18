#ifndef QUICKSORT_UTILS_H
#define QUICKSORT_UTILS_H

#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

void swap(int *e1, int *e2);
void merge_ascending(int *v1, int n1, int *v2, int n2, int *result);
void merge_descending(int *v1, int n1, int *v2, int n2, int *result);
void select_pivot_and_partition(int **elements_ptr, int n, MPI_Comm communicator,
                               int *num_small_local, int *num_large_local);


void determine_process_groups(int my_rank, int num_procs, int num_small_local, int num_large_local,
                              int **elements_ptr, int *group, int *partner_rank,
                              int **send_buf_ptr, int *send_count,
                              int **keep_buf_start_ptr, int *keep_count);

int* exchange_data(int *send_buf_ptr, int send_count, int partner_rank,
    MPI_Comm communicator, int my_rank, int *recv_count_ptr);   

int* merge_data(int *keep_buf_start_ptr, int keep_count, int *recv_buf, int recv_count,
                int my_rank, MPI_Comm communicator, int *new_n_ptr);

                
int recursive_sort(int **elements_ptr, int n, MPI_Comm communicator, int group, int my_rank, int pivot_strategy);

int global_sort(int **elements_ptr, int n, MPI_Comm communicator, int is_descending);


#endif /* QUICKSORT_UTILS_H */ 