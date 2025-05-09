// quicksort.c
#include "quicksort.h"
#include "pivot.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

int check_and_print(int *elements, int n, char *file_name) {
    FILE *fp = fopen(file_name, "w");
    if (!fp) return -2;
    
    for (int i = 0; i < n-1; i++) {
        if (elements[i] > elements[i+1]) {
            printf("Error: Elements %d and %d are out of order\n", elements[i], elements[i+1]);
            break;
        }
    }
    
    for (int i = 0; i < n; i++) fprintf(fp, "%d\n", elements[i]);
    fclose(fp);
    return 0;
}

int distribute_from_root(int *all_elements, int n, int **my_elements) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int *counts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));
    int remainder = n % size;
    
    for (int i = 0; i < size; i++) {
        counts[i] = n / size + (i < remainder);
        displs[i] = (i == 0) ? 0 : displs[i-1] + counts[i-1];
    }
    
    int local_count = counts[rank]; // 保存本地数据量
    *my_elements = malloc(local_count * sizeof(int));
    
    MPI_Scatterv(all_elements, counts, displs, MPI_INT,
                *my_elements, local_count, MPI_INT, ROOT, MPI_COMM_WORLD);
    
    free(counts);
    free(displs);
    
    return local_count; // 返回保存的值
}

void gather_on_root(int *all_elements, int *my_elements, int local_n) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int *counts = NULL, *displs = NULL;
    if (rank == ROOT) {
        counts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        MPI_Gather(&local_n, 1, MPI_INT, counts, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        
        displs[0] = 0;
        for (int i = 1; i < size; i++) 
            displs[i] = displs[i-1] + counts[i-1];
    } else {
        MPI_Gather(&local_n, 1, MPI_INT, NULL, 0, MPI_INT, ROOT, MPI_COMM_WORLD);
    }
    
    MPI_Gatherv(my_elements, local_n, MPI_INT, all_elements, counts, displs, 
               MPI_INT, ROOT, MPI_COMM_WORLD);
    
    if (rank == ROOT) {
        free(counts);
        free(displs);
    }
}

int global_sort(int **elements, int n, MPI_Comm comm, int pivot_strategy) {
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    
    if (size == 1) return n;
    
    // Select pivot and split data
    int split_idx = select_pivot(pivot_strategy, *elements, n, comm);
    // int pivot = (split_idx < n) ? (*elements)[split_idx] : (*elements)[n-1];
    
    // Split communicator
    int color = (rank < size/2) ? 0 : 1;
    MPI_Comm new_comm;
    MPI_Comm_split(comm, color, rank, &new_comm);
    
    // Exchange data between paired processes
    int partner_rank = rank + (color ? -size/2 : size/2);
    int send_size = color ? split_idx : (n - split_idx); // 交换条件
    int recv_size;
    
    MPI_Sendrecv(&send_size, 1, MPI_INT, partner_rank, 0,
                &recv_size, 1, MPI_INT, partner_rank, 0, comm, MPI_STATUS_IGNORE);
    
    int *send_buf = color ? *elements : (*elements + split_idx); // 修正缓冲区指针
    int *recv_buf = malloc(recv_size * sizeof(int));
    
    MPI_Sendrecv(send_buf, send_size, MPI_INT, partner_rank, 0,
                recv_buf, recv_size, MPI_INT, partner_rank, 0, comm, MPI_STATUS_IGNORE);
    
    // Merge and sort
    int new_n = (color ? (n - split_idx) : split_idx) + recv_size; // 修正本地保留长度
    int *merged = malloc(new_n * sizeof(int));

    if (color) {
        merge_ascending(*elements + split_idx, n - split_idx, recv_buf, recv_size, merged);
    } else {
        merge_ascending(*elements, split_idx, recv_buf, recv_size, merged);
    }
    
    free(*elements);
    free(recv_buf);
    *elements = merged;
    
    // Recursive call
    int final_n = global_sort(elements, new_n, new_comm, pivot_strategy);
    MPI_Comm_free(&new_comm);
    return final_n;
}

void merge_ascending(int *v1, int n1, int *v2, int n2, int *result) {
    int i = 0, j = 0, k = 0;
    while (i < n1 && j < n2) {
        if (v1[i] < v2[j]) result[k++] = v1[i++];
        else result[k++] = v2[j++];
    }
    while (i < n1) result[k++] = v1[i++];
    while (j < n2) result[k++] = v2[j++];
}

int read_input(char *file_name, int **elements) {
    FILE *fp = fopen(file_name, "r");
    if (!fp) return -1;
    
    int n, val, capacity = 100;
    *elements = malloc(capacity * sizeof(int));
    n = 0;
    
    while (fscanf(fp, "%d", &val) != EOF) {
        if (n >= capacity) {
            capacity *= 2;
            *elements = realloc(*elements, capacity * sizeof(int));
        }
        (*elements)[n++] = val;
    }
    fclose(fp);
    return n;
}

int sorted_ascending(int *elements, int n) {
    for (int i = 0; i < n-1; i++) {
        if (elements[i] > elements[i+1]) return 0;
    }
    return 1;
}

void swap(int *e1, int *e2) {
    int tmp = *e1;
    *e1 = *e2;
    *e2 = tmp;
}