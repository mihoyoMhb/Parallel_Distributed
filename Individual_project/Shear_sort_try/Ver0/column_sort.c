#include "column_sort.h"
#include <stdio.h>
#include <stdlib.h> // For qsort, malloc, free
#include <string.h> // For memcpy
#include <mpi.h>
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
 * Perform the column sorting phase of shearsort.
 * Columns are always sorted in ascending order.
 * If columns are distributed (p_r > 1), each global column is gathered across the relevant process column,
 * sorted, and then the appropriate segment is taken back by each process.
 * If columns are not distributed vertically (p_r == 1), each process sorts its local column segments directly
 * (though this data is non-contiguous and requires temporary buffering).
 */
int perform_column_sort(DataDistribution *distrib, int phase_idx, int pivot_strategy) {
    (void)pivot_strategy; // Not used with qsort-based approach
    (void)phase_idx;      // Not used

    if (!distrib || !distrib->local_data) {
        fprintf(stderr, "Process %d: Invalid DataDistribution in perform_column_sort.\n", distrib ? distrib->world_rank : -1);
        return -1;
    }

    int local_rows = distrib->local_rows;
    int local_cols = distrib->local_cols;
    int N = distrib->N;
    int my_rank_in_col_comm = -1;
    int procs_in_col_comm = -1; // Will be p_r

    if (N == 0 || local_rows == 0) return 0; // No data to sort in columns / no rows

    if (distrib->col_comm != MPI_COMM_NULL) {
        MPI_Comm_rank(distrib->col_comm, &my_rank_in_col_comm);
        MPI_Comm_size(distrib->col_comm, &procs_in_col_comm);
    }

    // Case 1: Columns are not split across process rows (each process holds full-height segments)
    // This happens if p_r = 1. Then procs_in_col_comm might be 1 or col_comm could be MPI_COMM_NULL.
    // A simpler check: if local_rows == N.
    if (local_rows == N) { // Or procs_in_col_comm <= 1
        int *temp_col_buffer = (int*)malloc(N * sizeof(int)); // N is local_rows here
        if (!temp_col_buffer) {
            fprintf(stderr, "Process %d: Failed to allocate temp_col_buffer in perform_column_sort (local_rows_N=%d)\n", distrib->world_rank, N);
            MPI_Abort(MPI_COMM_WORLD, 1); return -1;
        }
        for (int lc = 0; lc < local_cols; lc++) {
            // Extract column into temp_col_buffer
            for (int lr = 0; lr < N; lr++) { // N is local_rows
                temp_col_buffer[lr] = distrib->local_data[lr * local_cols + lc];
            }
            // Sort it
            qsort(temp_col_buffer, N, sizeof(int), compare_ascending_col);
            // Copy back
            for (int lr = 0; lr < N; lr++) {
                distrib->local_data[lr * local_cols + lc] = temp_col_buffer[lr];
            }
        }
        free(temp_col_buffer);
        return 0;
    }

    // Case 2: Columns are split across process rows (procs_in_col_comm > 1)
    if (procs_in_col_comm <= 0) {
        if (distrib->world_rank == 0) fprintf(stderr, "Process %d: Column communicator invalid or has no processes (%d) when expected for distributed column sort.\n", distrib->world_rank, procs_in_col_comm);
        MPI_Abort(MPI_COMM_WORLD, 1); return -1;
    }

    int *full_col_data = (int*)malloc(N * sizeof(int));
    if (!full_col_data) {
        fprintf(stderr, "Process %d: Failed to allocate memory for full_col_data (N=%d) in perform_column_sort\n", distrib->world_rank, N);
        MPI_Abort(MPI_COMM_WORLD, 1); return -1;
    }

    int *my_col_segment_send_buf = (int*)malloc(local_rows * sizeof(int));
    if (!my_col_segment_send_buf) {
        fprintf(stderr, "Process %d: Failed to allocate memory for my_col_segment_send_buf (local_rows=%d) in perform_column_sort\n", distrib->world_rank, local_rows);
        free(full_col_data);
        MPI_Abort(MPI_COMM_WORLD, 1); return -1;
    }

    for (int lc = 0; lc < local_cols; lc++) { // Iterate over local columns this process has a part of
        // Extract this process's segment of the current global column
        for (int lr = 0; lr < local_rows; lr++) {
            my_col_segment_send_buf[lr] = distrib->local_data[lr * local_cols + lc];
        }

        // Gather all segments of the current global column from processes in col_comm
        MPI_Allgather(my_col_segment_send_buf, local_rows, MPI_INT,
                      full_col_data, local_rows, MPI_INT, distrib->col_comm);
        
        // Sort the entire gathered global column (always ascending for column sort)
        qsort(full_col_data, N, sizeof(int), compare_ascending_col);

        // Copy the locally relevant sorted segment back to distrib->local_data
        int start_offset_in_full_col = my_rank_in_col_comm * local_rows;
        if (start_offset_in_full_col + local_rows <= N) {
            for (int lr = 0; lr < local_rows; lr++) {
                distrib->local_data[lr * local_cols + lc] = full_col_data[start_offset_in_full_col + lr];
            }
        } else {
             fprintf(stderr, "Process %d (Rank in ColComm %d, GridCol %d): Error - Column data copy bounds. Offset %d + local_rows %d > N %d for local_col_idx %d\n",
                distrib->world_rank, my_rank_in_col_comm, distrib->my_grid_col, start_offset_in_full_col, local_rows, N, lc);
            free(my_col_segment_send_buf);
            free(full_col_data);
            MPI_Abort(MPI_COMM_WORLD, 1); return -1;
        }
    }

    free(my_col_segment_send_buf);
    free(full_col_data);

    return 0;
}