#include "column_sort_v2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

// Comparison function for qsort (ascending)
static int compare_ascending_col(const void *a, const void *b) {
    int int_a = *((const int*)a);
    int int_b = *((const int*)b);
    if (int_a < int_b) return -1;
    if (int_a > int_b) return 1;
    return 0;
}

/**
 * Perform the column sorting phase of shearsort (V2).
 * Columns are always sorted in ascending order.
 * This version expects data_distribution_v2 (1D row decomposition, p_c=1).
 * Each process owns `local_rows` full-width segments of global columns.
 * It gathers these segments using MPI_Allgather, sorts the full column, then takes back its segment.
 */
int perform_column_sort_v2(DataDistribution *distrib, int phase_idx, int pivot_strategy) {
    (void)phase_idx;      // Not used in this qsort-based column sort
    (void)pivot_strategy; // Not used

    if (!distrib) {
        fprintf(stderr, "Process (unknown rank): NULL DataDistribution in perform_column_sort_v2.\\n");
        return -1;
    }
    if (!distrib->local_data && (distrib->local_rows * distrib->local_cols > 0) ){
        // Data expected but not present
        fprintf(stderr, "Process %d: NULL local_data in perform_column_sort_v2 when data expected (local_rows=%d, local_cols=%d).\\n", 
                distrib->world_rank, distrib->local_rows, distrib->local_cols);
        return -1;
    }
    // If local_rows or local_cols is 0, local_data might be NULL, which is fine.

    int local_rows = distrib->local_rows; // Number of rows this process owns a segment of for each column
    int local_cols = distrib->local_cols; // Should be N (total number of columns in the matrix)
    int N = distrib->N;                   // Global matrix dimension (N x N)

    if (N == 0) return 0; // No data to sort in a 0x0 matrix
    if (local_rows == 0) return 0; // This process has no rows, so no column segments to sort
    if (local_cols == 0) return 0; // Should not happen if N > 0, as local_cols = N

    // With 1D row decomposition (p_c=1):
    // distrib->col_comm should include all processes.
    // distrib->my_grid_col is 0 for all processes.
    // distrib->p_r is world_size.
    // local_cols is N.

    int my_rank_in_col_comm = -1;
    int procs_in_col_comm = -1;

    if (distrib->col_comm == MPI_COMM_NULL) {
        if (distrib->world_size > 1) { // If more than one proc, col_comm should be valid
            fprintf(stderr, "Process %d: Column communicator is MPI_COMM_NULL unexpectedly (world_size=%d).\\n", 
                    distrib->world_rank, distrib->world_size);
            return -1; 
        } else { // Single process case
            // Fall through to local sort logic (adapted from original Case 1)
             my_rank_in_col_comm = 0;
             procs_in_col_comm = 1;
        }
    } else {
        MPI_Comm_rank(distrib->col_comm, &my_rank_in_col_comm);
        MPI_Comm_size(distrib->col_comm, &procs_in_col_comm);
    }

    // The original "Case 1" (local_rows == N) implies a single process owns the entire column height.
    // In our 1D row decomposition (p_c=1), local_rows = N / world_size.
    // So, local_rows == N only if world_size == 1.
    if (procs_in_col_comm == 1) { // Equivalent to world_size == 1 for our 1D setup
        if (local_rows != N) { // Sanity check for single process case
             fprintf(stderr, "Process %d: In single proc mode for column sort, but local_rows (%d) != N (%d).\\n", 
                     distrib->world_rank, local_rows, N);
             return -1; // Should not happen
        }
        int *temp_col_buffer = (int*)malloc(N * sizeof(int)); // N is local_rows here
        if (!temp_col_buffer) {
            fprintf(stderr, "Process %d: Failed to allocate temp_col_buffer (N=%d) for local column sort.\\n", distrib->world_rank, N);
            MPI_Abort(MPI_COMM_WORLD, 1); return -1;
        }
        for (int lc = 0; lc < local_cols; lc++) { // local_cols is N here
            for (int lr = 0; lr < N; lr++) { // lr iterates through the full column height (local_rows)
                temp_col_buffer[lr] = distrib->local_data[lr * local_cols + lc];
            }
            qsort(temp_col_buffer, N, sizeof(int), compare_ascending_col);
            for (int lr = 0; lr < N; lr++) {
                distrib->local_data[lr * local_cols + lc] = temp_col_buffer[lr];
            }
        }
        free(temp_col_buffer);
        return 0;
    }

    // Case 2: Columns are split across processes (procs_in_col_comm > 1, i.e. world_size > 1)
    // This is the expected path for multi-process runs with 1D row decomposition.
    if (procs_in_col_comm <= 0) { // Should have been caught if col_comm was NULL and world_size > 1
        fprintf(stderr, "Process %d: Column communicator invalid or has %d processes when >1 expected.\\n", 
                distrib->world_rank, procs_in_col_comm);
        MPI_Abort(MPI_COMM_WORLD, 1); return -1;
    }
    if (distrib->col_comm == MPI_COMM_NULL && procs_in_col_comm > 1) {
         fprintf(stderr, "Process %d: col_comm is NULL but procs_in_col_comm = %d. Contradiction.\\n", distrib->world_rank, procs_in_col_comm);
         MPI_Abort(MPI_COMM_WORLD, 1); return -1;
    }


    int *full_col_data_buffer = (int*)malloc(N * sizeof(int)); // Buffer for one full column
    if (!full_col_data_buffer) {
        fprintf(stderr, "Process %d: Failed to allocate full_col_data_buffer (N=%d) in perform_column_sort_v2.\\n", distrib->world_rank, N);
        MPI_Abort(MPI_COMM_WORLD, 1); return -1;
    }

    // Each process has `local_rows` elements for each of the `local_cols` (which is N) columns.
    // This send buffer is for one column segment from this process.
    int *my_col_segment_send_buf = (int*)malloc(local_rows * sizeof(int)); 
    if (!my_col_segment_send_buf) {
        fprintf(stderr, "Process %d: Failed to allocate my_col_segment_send_buf (local_rows=%d).\\n", distrib->world_rank, local_rows);
        free(full_col_data_buffer);
        MPI_Abort(MPI_COMM_WORLD, 1); return -1;
    }

    // Iterate over all columns of the matrix (local_cols is N)
    for (int global_col_idx = 0; global_col_idx < N; global_col_idx++) {
        // Extract this process's segment of the current global_col_idx
        for (int lr = 0; lr < local_rows; lr++) {
            // local_data is [row][col], so lr * N + global_col_idx
            my_col_segment_send_buf[lr] = distrib->local_data[lr * N + global_col_idx];
        }

        // Gather all segments of the current global column using the column communicator
        // Each of `procs_in_col_comm` (which is world_size) processes sends `local_rows` elements.
        // The result `full_col_data_buffer` will have `procs_in_col_comm * local_rows = world_size * (N/world_size) = N` elements.
        MPI_Allgather(my_col_segment_send_buf, local_rows, MPI_INT,
                      full_col_data_buffer, local_rows, MPI_INT, distrib->col_comm);
        
        // Sort the entire gathered global column (always ascending for column sort)
        qsort(full_col_data_buffer, N, sizeof(int), compare_ascending_col);

        // Copy the locally relevant sorted segment back to distrib->local_data
        // my_rank_in_col_comm is this process's rank within the group of processes sharing this column (i.e. its world_rank)
        // Its segment starts at offset `my_rank_in_col_comm * local_rows` in the `full_col_data_buffer`.
        int start_offset_in_full_col = my_rank_in_col_comm * local_rows;
        
        // Boundary check, should always pass if N is multiple of world_size and logic is correct
        if (start_offset_in_full_col + local_rows > N) {
             fprintf(stderr, "Process %d (Rank in ColComm %d): Error - Column data copy bounds. Offset %d + local_rows %d > N %d for global_col_idx %d\\n",
                distrib->world_rank, my_rank_in_col_comm, start_offset_in_full_col, local_rows, N, global_col_idx);
            free(my_col_segment_send_buf);
            free(full_col_data_buffer);
            MPI_Abort(MPI_COMM_WORLD, 1); return -1;
        }

        for (int lr = 0; lr < local_rows; lr++) {
            distrib->local_data[lr * N + global_col_idx] = full_col_data_buffer[start_offset_in_full_col + lr];
        }
    }

    free(my_col_segment_send_buf);
    free(full_col_data_buffer);

    return 0;
} 