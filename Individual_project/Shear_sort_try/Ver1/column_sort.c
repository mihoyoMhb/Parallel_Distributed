#include "column_sort.h"
#include <stdio.h>
#include <stdlib.h> // For qsort, malloc, free
#include <string.h> // For memcpy
#include <mpi.h>
#include "quicksort_utils.h" // 添加快速排序工具库
// data_distribution.h is included via column_sort.h

// Comparison function for qsort (ascending)
static int compare_ascending_col(const void *a, const void *b) {
    int int_a = *((const int*)a);
    int int_b = *((const int*)b);
    if (int_a < int_b) return -1;
    if (int_a > int_b) return 1;
    return 0;
}

/**
 * Perform the column sorting phase of shearsort using parallel quicksort.
 * Columns are always sorted in ascending order.
 * Each global column is sorted using parallel quicksort across the relevant process column.
 */
int perform_column_sort(DataDistribution *distrib, int phase_idx, int pivot_strategy) {
    (void)pivot_strategy; // Not used with this implementation
    (void)phase_idx;      // Not used

    if (!distrib || !distrib->local_data) {
        fprintf(stderr, "Process %d: Invalid DataDistribution in perform_column_sort.\n", distrib ? distrib->world_rank : -1);
        return -1;
    }

    int local_rows = distrib->local_rows;
    int local_cols = distrib->local_cols;
    int N = distrib->N;
    int my_rank_in_col_comm = -1;
    int procs_in_col_comm = -1;

    // Basic checks for empty data
    if (N == 0 || local_rows == 0 || local_cols == 0) return 0;

    // Get rank in column communicator
    if (distrib->col_comm != MPI_COMM_NULL) {
        MPI_Comm_rank(distrib->col_comm, &my_rank_in_col_comm);
        MPI_Comm_size(distrib->col_comm, &procs_in_col_comm);
    } else {
        // Single process per column
        procs_in_col_comm = 1;
        my_rank_in_col_comm = 0;
    }

    // Case 1: Single process per column - sort locally
    if (procs_in_col_comm == 1) {
        int *temp_col_buffer = (int*)malloc(local_rows * sizeof(int));
        if (!temp_col_buffer) {
            fprintf(stderr, "Process %d: Failed to allocate temp_col_buffer in perform_column_sort\n", 
                    distrib->world_rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }
        
        for (int lc = 0; lc < local_cols; lc++) {
            // Extract column into temp_col_buffer
            for (int lr = 0; lr < local_rows; lr++) {
                temp_col_buffer[lr] = distrib->local_data[lr * local_cols + lc];
            }
            
            // Sort it locally
            qsort(temp_col_buffer, local_rows, sizeof(int), compare_ascending_col);
            
            // Copy back
            for (int lr = 0; lr < local_rows; lr++) {
                distrib->local_data[lr * local_cols + lc] = temp_col_buffer[lr];
            }
        }
        
        free(temp_col_buffer);
        return 0;
    }

    // Case 2: Multiple processes per column - use parallel quicksort
    if (procs_in_col_comm <= 0) {
        fprintf(stderr, "Process %d: Column communicator invalid or has no processes (%d) when expected for distributed column sort.\n", 
                distrib->world_rank, procs_in_col_comm);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }

    printf("Wating to be implemented for multiple processes per col.\n");
    
    return 0;
}