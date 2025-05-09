#include "quicksort.h"
#include "pivot.h" 
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For memcpy

// Swap two integers
void swap(int *e1, int *e2) {
    int temp = *e1;
    *e1 = *e2;
    *e2 = temp;
}

// Check if elements are sorted in ascending order
int sorted_ascending(int *elements, int n) {
    if (n == 0 || elements == NULL) return 1; // Empty or NULL array is considered sorted
    for (int i = 0; i < n - 1; ++i) {
        if (elements[i] > elements[i + 1]) {
            // fprintf(stderr, "Error: Array not sorted. elements[%d]=%d > elements[%d]=%d\n",
            //         i, elements[i], i + 1, elements[i + 1]);
            return 0; // Not sorted
        }
    }
    return 1; // Sorted
}

// Merge two sorted arrays v1 and v2 into result
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

// Read problem size and elements from file
int read_input(char *file_name, int **elements_ptr) {
    FILE *file = fopen(file_name, "r");
    if (!file) {
        perror("Failed to open input file");
        *elements_ptr = NULL;
        return -1; 
    }

    int n;
    if (fscanf(file, "%d", &n) != 1) {
        fprintf(stderr, "Failed to read problem size from file %s\n", file_name);
        fclose(file);
        *elements_ptr = NULL;
        return -1; 
    }

    if (n < 0) {
        fprintf(stderr, "Problem size n cannot be negative (%d)\n", n);
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
    if (!*elements_ptr) {
        fprintf(stderr, "Failed to allocate memory for %d elements in read_input\n", n);
        fclose(file);
        return -1; 
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

// Verify sorted and print to file
int check_and_print(int *elements, int n, char *file_name) {
    if (!sorted_ascending(elements, n)) {
        fprintf(stderr, "Warning: Final array is not sorted. Printing to file %s anyway.\n", file_name);
    }

    FILE *file = fopen(file_name, "w");
    if (!file) {
        perror("Failed to open output file for writing");
        return -2; // I/O error as per spec
    }

    for (int i = 0; i < n; ++i) {
        fprintf(file, "%d", elements[i]);
        if (i < n - 1) {
            fprintf(file, " "); 
        }
    }
    fprintf(file, "\n"); 

    fclose(file);
    return 0; // Success
}

// Distribute elements from root to all processes
int distribute_from_root(int *all_elements_on_root, int n_total, int **my_elements_ptr) {
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int *sendcounts = NULL;
    int *displs = NULL;
    int local_n;

    if (rank == ROOT) {
        if (comm_size > 0) {
            sendcounts = (int *)calloc(comm_size, sizeof(int)); // Use calloc for safety
            displs = (int *)calloc(comm_size, sizeof(int));
            if (!sendcounts || !displs) {
                 fprintf(stderr, "ROOT: Failed to allocate memory for scatterv metadata.\n");
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
    
    if (local_n < 0) { 
        fprintf(stderr, "[Process %d] Received negative local_n = %d in distribute_from_root.\n", rank, local_n);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (local_n > 0) {
        *my_elements_ptr = (int *)malloc(local_n * sizeof(int));
        if (!*my_elements_ptr) {
            fprintf(stderr, "[Process %d] Failed to allocate memory for %d local elements.\n", rank, local_n);
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

// Gather elements from all processes on root
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
                 MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    MPI_Gather(&local_n, 1, MPI_INT, recvcounts, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
        if (comm_size > 0) {
            int current_displ = 0;
            for (int i = 0; i < comm_size; ++i) {
                displs[i] = current_displ;
                current_displ += recvcounts[i];
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

// Perform the global part of parallel quick sort
int global_sort(int **elements_ptr, int n, MPI_Comm communicator, int pivot_strategy) {
    int my_rank, num_procs;
    MPI_Comm_rank(communicator, &my_rank);
    MPI_Comm_size(communicator, &num_procs);

    if (num_procs <= 1) {
        return n; 
    }

    // Step 3.1 & 3.2 (implicit): Select pivot, get local split index k
    int k = select_pivot(pivot_strategy, *elements_ptr, n, communicator);
    int num_small_local = k;
    int num_large_local = n - k;

    // Step 3.3: Split processes and exchange data
    // This assumes num_procs is a power of 2 for clean halving.
    int color = (my_rank < num_procs / 2) ? 0 : 1; 
    int partner_rank;
    int *send_buf_ptr;
    int send_count;
    int *keep_buf_start_ptr; 
    int keep_count;

    if (color == 0) { // Lower group (ranks 0 to num_procs/2 - 1)
        partner_rank = my_rank + num_procs / 2;
        send_buf_ptr = (*elements_ptr) + num_small_local; 
        send_count = num_large_local;
        keep_buf_start_ptr = *elements_ptr; 
        keep_count = num_small_local;
    } else { // Upper group (ranks num_procs/2 to num_procs - 1)
        partner_rank = my_rank - num_procs / 2;
        send_buf_ptr = *elements_ptr; 
        send_count = num_small_local;
        keep_buf_start_ptr = (*elements_ptr) + num_small_local; 
        keep_count = num_large_local;
    }

    int recv_count;
    MPI_Sendrecv(&send_count, 1, MPI_INT, partner_rank, 0, // Tag 0 for counts
                 &recv_count, 1, MPI_INT, partner_rank, 0,
                 communicator, MPI_STATUS_IGNORE);

    int *recv_buf = NULL;
    if (recv_count > 0) {
        recv_buf = (int *)malloc(recv_count * sizeof(int));
        if (!recv_buf) {
            fprintf(stderr, "[Proc %d] Failed to alloc recv_buf for %d ints\n", my_rank, recv_count);
            MPI_Abort(communicator, 1);
        }
    }
    
    MPI_Sendrecv(send_buf_ptr, send_count, MPI_INT, partner_rank, 1, // Tag 1 for data
                 recv_buf, recv_count, MPI_INT, partner_rank, 1,
                 communicator, MPI_STATUS_IGNORE);

    // Step 3.4: Merge kept and received data
    int new_n = keep_count + recv_count;
    int *new_elements_array = NULL;

    if (new_n > 0) {
        new_elements_array = (int *)malloc(new_n * sizeof(int));
        if (!new_elements_array) {
            fprintf(stderr, "[Proc %d] Failed to alloc new_elements_array for %d ints\n", my_rank, new_n);
            MPI_Abort(communicator, 1);
        }

        if (keep_count > 0 && recv_count > 0) {
            // Need a distinct copy of the kept part if merging in-place is tricky or keep_buf_start_ptr is part of *elements_ptr
            // For simplicity, copy kept part to a temp buffer for merging if it's not empty.
            int* temp_kept_elements = (int*)malloc(keep_count * sizeof(int));
             if (!temp_kept_elements) {
                 fprintf(stderr, "[Proc %d] Failed to alloc temp_kept_elements for %d ints\n", my_rank, keep_count);
                 MPI_Abort(communicator, 1);
            }
            memcpy(temp_kept_elements, keep_buf_start_ptr, keep_count * sizeof(int));
            merge_ascending(temp_kept_elements, keep_count, recv_buf, recv_count, new_elements_array);
            free(temp_kept_elements);
        } else if (keep_count > 0) { 
            memcpy(new_elements_array, keep_buf_start_ptr, keep_count * sizeof(int));
        } else if (recv_count > 0) { 
            memcpy(new_elements_array, recv_buf, recv_count * sizeof(int));
        }
        // If both keep_count and recv_count are 0, new_n is 0, this block is skipped.
    }
    
    if (*elements_ptr) { 
        free(*elements_ptr);
    }
    *elements_ptr = new_elements_array; // Can be NULL if new_n is 0
    n = new_n;

    if (recv_buf) {
        free(recv_buf);
    }

    // Step 4: Repeat recursively for each half
    MPI_Comm new_comm;
    MPI_Comm_split(communicator, color, my_rank, &new_comm);
    n = global_sort(elements_ptr, n, new_comm, pivot_strategy);
    MPI_Comm_free(&new_comm);
    
    return n;
}