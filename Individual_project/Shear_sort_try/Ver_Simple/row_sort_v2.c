#include "row_sort_v2.h"
#include <stdio.h>
#include <stdlib.h> // For qsort
#include <mpi.h>    // For MPI_Barrier (optional, for debugging or strict phase separation)

// Comparison function for qsort (ascending)
static int compare_ascending(const void *a, const void *b) {
    int int_a = *((const int*)a);
    int int_b = *((const int*)b);
    if (int_a < int_b) return -1;
    if (int_a > int_b) return 1;
    return 0;
}

// Comparison function for qsort (descending)
static int compare_descending(const void *a, const void *b) {
    int int_a = *((const int*)a);
    int int_b = *((const int*)b);
    if (int_a > int_b) return -1; // Note: For descending, larger value comes first
    if (int_a < int_b) return 1;
    return 0;
}

int perform_row_sort_v2(DataDistribution *distrib, int phase_idx) {
    if (!distrib || !distrib->local_data) {
        if (distrib && distrib->local_rows * distrib->local_cols > 0) { // Only error if data was expected
            fprintf(stderr, "Process %d: Invalid DataDistribution or NULL local_data in perform_row_sort_v2.\\n", 
                    distrib ? distrib->world_rank : -1);
            return -1;
        }
        if (!distrib) return -1; // distrib itself is NULL
        // If local_data is NULL but block size is 0, it's not an error for sorting.
    }

    if (distrib->local_rows == 0 || distrib->local_cols == 0) {
        return 0; // No data to sort for this process
    }

    // For 1D row decomposition, distrib->my_grid_row is essentially the process rank.
    // distrib->local_rows is the number of full rows this process owns.
    // distrib->local_cols is N (the full width of the matrix).

    int start_global_row_for_process = distrib->my_grid_row * distrib->local_rows;

    for (int lr = 0; lr < distrib->local_rows; lr++) {
        int current_global_row_idx = start_global_row_for_process + lr;
        int *row_data_ptr = &distrib->local_data[lr * distrib->local_cols];

        // Determine sort direction based on global row index (0-indexed)
        // Even global rows (0, 2, ...) -> ascending
        // Odd global rows (1, 3, ...) -> descending
        if (current_global_row_idx % 2 == 0) {
            qsort(row_data_ptr, distrib->local_cols, sizeof(int), compare_ascending);
        } else {
            qsort(row_data_ptr, distrib->local_cols, sizeof(int), compare_descending);
        }
    }
    
    // Optional: Barrier to ensure all row sorts complete before proceeding (e.g., to column sort)
    // Useful for debugging or strict phase synchronization, though Shear Sort already has barriers between phases.
    // MPI_Barrier(MPI_COMM_WORLD); 

    return 0;
} 