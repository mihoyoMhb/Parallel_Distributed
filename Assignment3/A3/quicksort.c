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
    // Assume file opens successfully in controlled tests

    int n;
    if (fscanf(file, "%d", &n) != 1) {
        fprintf(stderr, "Error reading!!!\n");
        fclose(file);
        *elements_ptr = NULL;
        return -1;
    }

    if (n == 0) {
        *elements_ptr = NULL;
        fclose(file);
        return 0;
    }

    *elements_ptr = malloc(n * sizeof(int));
    if (!*elements_ptr) {
        fprintf(stderr, "Error reading!!!\n");
        fclose(file);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < n; ++i) {
        if (fscanf(file, "%d", &((*elements_ptr)[i])) != 1) {
            fprintf(stderr, "Error reading!!! %d\n", i);
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
    // Assuming sorted_ascending has no side effects besides checking
    sorted_ascending(elements, n);

    FILE *file = fopen(file_name, "w");

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

    int *sendcounts = NULL, *displacements = NULL, local_n;

    if (rank == ROOT) {
        sendcounts = calloc(comm_size, sizeof(int));
        displacements     = calloc(comm_size, sizeof(int));
        // Still the same as previous assignment 1 2
        // Distribute elements evenly across processes
        int base    = n_total / comm_size;
        int remainder     = n_total % comm_size;
        int offset  = 0;
        for (int i = 0; i < comm_size; i++) {
            if (i < remainder) {
                sendcounts[i] = base + 1;
            } else {
                sendcounts[i] = base;
            }
            displacements[i]     = offset;
            offset       += sendcounts[i];
        }
    }

    MPI_Scatter(sendcounts, 1, MPI_INT,
                &local_n,    1, MPI_INT,
                ROOT, MPI_COMM_WORLD);
    *my_elements_ptr = malloc(local_n * sizeof(int));
    // Scatter the elements to all processes
    MPI_Scatterv(all_elements_on_root, sendcounts, displacements, MPI_INT,
                 *my_elements_ptr, local_n, MPI_INT,
                 ROOT, MPI_COMM_WORLD);

    free(sendcounts);
    free(displacements);
    return local_n;
}

void gather_on_root(int *all_elements_buffer_on_root, int *my_elements, int local_n) {
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int *recvcounts = NULL, *displacements = NULL;

    if (rank == ROOT) {
        recvcounts = malloc(comm_size * sizeof(int));
        displacements     = malloc(comm_size * sizeof(int));
    }

    MPI_Gather(&local_n, 1, MPI_INT,
               recvcounts, 1, MPI_INT,
               ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
        int offset = 0;
        for (int i = 0; i < comm_size; i++) {
            displacements[i] = offset;
            offset   += recvcounts[i];
        }
    }

    MPI_Gatherv(my_elements, local_n, MPI_INT,
                all_elements_buffer_on_root, recvcounts, displacements, MPI_INT,
                ROOT, MPI_COMM_WORLD);

    free(recvcounts);
    free(displacements);
}

/**
 * Select pivot and partition locally.
 * 
 * @param elements_ptr    pointer to the elements array
 * @param n               number of elements
 * @param communicator    MPI communicator
 * @param pivot_strategy  pivot selection strategy
 * @param num_small_local returns count of <= pivot elements
 * @param num_large_local returns count of > pivot elements
 * @return number of <= pivot elements
 */
int select_pivot_and_partition(int **elements_ptr, int n, MPI_Comm communicator,
                               int pivot_strategy, int *num_small_local, int *num_large_local) {
    int k = select_pivot(pivot_strategy, *elements_ptr, n, communicator);
    *num_small_local = k;
    *num_large_local = n - k;
    return k;
}

/**
 * Determine process group, partner, and send/keep buffers.
 * 
 * @param my_rank            rank of this process
 * @param num_procs          total number of processes
 * @param num_small_local    count of <= pivot
 * @param num_large_local    count of > pivot
 * @param elements_ptr       pointer to elements array
 * @param group              returns group color (0=low,1=high)
 * @param partner_rank       returns partner rank
 * @param send_buf_ptr       returns pointer to send buffer start
 * @param send_count         returns number of elements to send
 * @param keep_buf_start_ptr returns pointer to keep buffer start
 * @param keep_count         returns number of elements to keep
 */
void determine_process_groups(int my_rank, int num_procs, int num_small_local, int num_large_local,
                              int **elements_ptr, int *group, int *partner_rank,
                              int **send_buf_ptr, int *send_count,
                              int **keep_buf_start_ptr, int *keep_count) {
    *group = (my_rank < num_procs / 2) ? 0 : 1;

    if (*group == 0) {
        *partner_rank      = my_rank + num_procs / 2;
        *send_buf_ptr      = (*elements_ptr) + num_small_local;
        *send_count        = num_large_local;
        *keep_buf_start_ptr= *elements_ptr;
        *keep_count        = num_small_local;
    } else {
        *partner_rank      = my_rank - num_procs / 2;
        *send_buf_ptr      = *elements_ptr;
        *send_count        = num_small_local;
        *keep_buf_start_ptr= (*elements_ptr) + num_small_local;
        *keep_count        = num_large_local;
    }
}

/**
 * Exchange data with partner.
 * 
 * @param send_buf_ptr  pointer to send buffer
 * @param send_count    number of elements to send
 * @param partner_rank  rank of partner process
 * @param communicator  MPI communicator
 * @param my_rank       rank of this process
 * @param recv_count_ptr returns number of elements received
 * @return pointer to received buffer (NULL if none)
 */
int* exchange_data(int *send_buf_ptr, int send_count, int partner_rank,
                   MPI_Comm communicator, int my_rank, int *recv_count_ptr) {
    int recv_count;
    // Exchage send_count with partner to determine recv_count
    MPI_Sendrecv(&send_count,  1, MPI_INT, partner_rank, 0,
                 &recv_count,  1, MPI_INT, partner_rank, 0,
                 communicator, MPI_STATUS_IGNORE);

    *recv_count_ptr = recv_count;
    int *recv_buf = NULL;
    
    // Allocate buffer for received data if recv_count > 0
    if (recv_count > 0) {
        recv_buf = malloc(recv_count * sizeof(int));
        if (!recv_buf) {
            MPI_Abort(communicator, 1);
        }
    }
    // Exchange actual data
    MPI_Sendrecv(send_buf_ptr, send_count, MPI_INT,  partner_rank, 1,
                 recv_buf,    recv_count, MPI_INT, partner_rank, 1,
                 communicator, MPI_STATUS_IGNORE);

    return recv_buf;
}

/**
 * Merge local kept data with received data.
 * 
 * @param keep_buf_start_ptr pointer to kept data buffer
 * @param keep_count         number of kept elements
 * @param recv_buf           received data buffer
 * @param recv_count         number of received elements
 * @param my_rank            rank of this process
 * @param communicator       MPI communicator
 * @param new_n_ptr          returns new total element count
 * @return pointer to merged array
 */
int* merge_data(int *keep_buf_start_ptr, int keep_count, int *recv_buf, int recv_count,
                int my_rank, MPI_Comm communicator, int *new_n_ptr) {
    int new_n = keep_count + recv_count;
    *new_n_ptr = new_n;
    int *new_elements_array = NULL;

    if (new_n > 0) {
        // Allocate new array for merged data
        new_elements_array = malloc(new_n * sizeof(int));
        if (keep_count > 0 && recv_count > 0) {
            // Merge two sorted arrays
            int *temp_array = malloc(keep_count * sizeof(int));
            if(!temp_array || !new_elements_array) {
                fprintf(stderr, "Memory allocation error in merge_data\n");
                MPI_Abort(communicator, 1);
            }
            memcpy(temp_array, keep_buf_start_ptr, keep_count * sizeof(int));
            merge_ascending(temp_array, keep_count, recv_buf, recv_count, new_elements_array);
            free(temp_array);
        } else if (keep_count > 0) {
            // Only kept data
            memcpy(new_elements_array, keep_buf_start_ptr, keep_count * sizeof(int));
        } else if (recv_count > 0) {
            // Only received data
            memcpy(new_elements_array, recv_buf, recv_count * sizeof(int));
        }
    }

    return new_elements_array;
}

/**
 * Split communicator and recurse.
 * 
 * @param elements_ptr    pointer to elements array
 * @param n               number of elements
 * @param communicator    current MPI communicator
 * @param group           group (0 or 1)
 * @param my_rank         rank of this process
 * @param pivot_strategy  pivot selection strategy
 * @return new element count after sort
 */
int recursive_sort(int **elements_ptr, int n, MPI_Comm communicator,
                   int group, int my_rank, int pivot_strategy) {
    // Create new communicator for the group
    MPI_Comm new_comm;
    MPI_Comm_split(communicator, group, my_rank, &new_comm);

    // Sort in the new communicator
    n = global_sort(elements_ptr, n, new_comm, pivot_strategy);
    // Free the new communicator
    MPI_Comm_free(&new_comm);
    return n;
}

/**
 * Main parallel quicksort.
 * 
 * @param elements_ptr    pointer to elements array
 * @param n               number of elements
 * @param communicator    MPI communicator
 * @param pivot_strategy  pivot selection strategy
 * @return sorted element count
 */
int global_sort(int **elements_ptr, int n, MPI_Comm communicator, int pivot_strategy) {
    int my_rank, num_procs;
    MPI_Comm_rank(communicator, &my_rank);
    MPI_Comm_size(communicator, &num_procs);

    // Base case: if only one process, no need to sort
    if (num_procs <= 1) {
        return n;
    }

    // 1. Select pivot and partition locally
    int num_small_local, num_large_local;
    select_pivot_and_partition(elements_ptr, n, communicator,
                               pivot_strategy, &num_small_local, &num_large_local);

    // 2. Determine groups and buffers
    int group, partner_rank;
    int *send_buf_ptr, *keep_buf_start_ptr;
    int send_count, keep_count;
    determine_process_groups(my_rank, num_procs, num_small_local, num_large_local,
                             elements_ptr, &group, &partner_rank,
                             &send_buf_ptr, &send_count,
                             &keep_buf_start_ptr, &keep_count);

    // 3. Exchange data
    int recv_count;
    int *recv_buf = exchange_data(send_buf_ptr, send_count, partner_rank,
                                  communicator, my_rank, &recv_count);

    // 4. Merge data
    int new_n;
    int *new_elements_array = merge_data(keep_buf_start_ptr, keep_count,
                                         recv_buf, recv_count,
                                         my_rank, communicator, &new_n);

    // Free old elements array and assign new one
    free(*elements_ptr);
    *elements_ptr = new_elements_array;
    n = new_n;

    // Free received buffer if it was allocated
    if (recv_buf) {
        free(recv_buf);
    }

    // 5. Recurse on new communicator
    n = recursive_sort(elements_ptr, n, communicator, group, my_rank, pivot_strategy);

    return n;
}
