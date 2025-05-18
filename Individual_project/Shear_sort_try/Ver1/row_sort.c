#include "row_sort.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "quicksort_utils.h"

#define ROOT 0

/**
 * Helper function to swap two integers
 */
void swap(int *e1, int *e2) {
    int temp = *e1;
    *e1 = *e2;
    *e2 = temp;
}

/**
 * Check if array is sorted in ascending order
 */
int sorted_ascending(int *elements, int n) {
    if (n == 0 || elements == NULL) return 1;
    for (int i = 0; i < n - 1; ++i) {
        if (elements[i] > elements[i + 1]) {
            return 0;
        }
    }
    return 1;
}

/**
 * Check if array is sorted in descending order
 */
int sorted_descending(int *elements, int n) {
    if (n == 0 || elements == NULL) return 1;
    for (int i = 0; i < n - 1; ++i) {
        if (elements[i] < elements[i + 1]) {
            return 0;
        }
    }
    return 1;
}

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
    // For descending, if a > b, a should come before b, so return -1
    if (int_a > int_b) return -1;
    if (int_a < int_b) return 1;
    return 0;
}

/**
 * Perform the row sorting phase of shearsort using parallel quicksort.
 * Each global row is gathered across the relevant process row,
 * sorted by parallel quicksort, and then the appropriate segment is taken back.
 */
int perform_row_sort(DataDistribution *distrib, int phase_idx, int pivot_strategy) {
    (void)pivot_strategy; // Not directly used in this implementation
    (void)phase_idx;      // Not directly used; sort direction based on global row index

    if (!distrib) {
        fprintf(stderr, "Process (unknown): Invalid DataDistribution (NULL) in perform_row_sort.\n");
        return -1;
    }
    if (!distrib->local_data) {
        fprintf(stderr, "Process %d: Invalid local_data (NULL) in perform_row_sort.\n", distrib->world_rank);
        return -1;
    }
    
    int local_rows = distrib->local_rows;
    int local_cols = distrib->local_cols;
    int N = distrib->N;
    int my_rank_in_row_comm = -1;
    int procs_in_row_comm = -1; 

    // Basic checks for empty data
    if (N == 0) return 0; 
    if (local_cols == 0 || local_rows == 0) return 0;

    // Get rank in row communicator
    if (distrib->row_comm != MPI_COMM_NULL) {
        MPI_Comm_rank(distrib->row_comm, &my_rank_in_row_comm);
        MPI_Comm_size(distrib->row_comm, &procs_in_row_comm);
    } else {
        // Single process per row
        procs_in_row_comm = 1; 
        my_rank_in_row_comm = 0; 
    }

    // Case 1: Single process per row - sort locally
    if (procs_in_row_comm == 1) { 
        for (int lr = 0; lr < local_rows; lr++) {
            int global_row_idx = distrib->my_grid_row * local_rows + lr;
            int is_descending = global_row_idx % 2; // Even rows asc, odd rows desc
            
            // If sorting in descending order, negate all values, sort, then negate back
            if (is_descending) {
                for (int j = 0; j < local_cols; j++) {
                    distrib->local_data[lr * local_cols + j] = -distrib->local_data[lr * local_cols + j];
                }
            }
            
            // Sort locally
            qsort(&(distrib->local_data[lr * local_cols]), local_cols, sizeof(int), compare_ascending);
            
            if (is_descending) {
                for (int j = 0; j < local_cols; j++) {
                    distrib->local_data[lr * local_cols + j] = -distrib->local_data[lr * local_cols + j];
                }
            }
        }
        return 0;
    }

    // Case 2: Multiple processes per row - use parallel quicksort
    if (procs_in_row_comm <= 0) { 
        fprintf(stderr, "Process %d: Row communicator has invalid size (%d) for distributed row sort.\n", 
                distrib->world_rank, procs_in_row_comm);
        MPI_Abort(MPI_COMM_WORLD, 1); 
        return -1;
    }

    // Process each row separately
    // for (int lr = 0; lr < local_rows; lr++) {
    //     int global_row_idx = distrib->my_grid_row * local_rows + lr;
    //     int is_descending = global_row_idx % 2; // Even rows asc, odd rows desc
        
    //     // Allocate a buffer for this row
    //     int *row_buffer = NULL;
    //     if (local_cols > 0) {
    //         row_buffer = (int*)malloc(local_cols * sizeof(int));
    //         if (!row_buffer) {
    //             fprintf(stderr, "Process %d: Failed to allocate memory for row_buffer in perform_row_sort\n", 
    //                     distrib->world_rank);
    //             MPI_Abort(MPI_COMM_WORLD, 1); 
    //             return -1;
    //         }
            
    //         // Copy local row segment to buffer
    //         for (int j = 0; j < local_cols; j++) {
    //             row_buffer[j] = distrib->local_data[lr * local_cols + j];
    //         }
            
    //         // If sorting in descending order, negate all values
    //         if (is_descending) {
    //             for (int j = 0; j < local_cols; j++) {
    //                 row_buffer[j] = -row_buffer[j];
    //             }
    //         }
    //     }
        
    //     // Sort the row using parallel quicksort (always use ascending)
    //     int row_size = local_cols;
    //     if (row_size > 0) {
    //         // Make sure all processes are synchronized before sort
    //         MPI_Barrier(distrib->row_comm);
            
    //         row_size = parallel_quicksort(&row_buffer, row_size, distrib->row_comm, 0);
            
    //         // If sorting in descending order, negate values back
    //         if (is_descending && row_size > 0 && row_buffer) {
    //             for (int j = 0; j < row_size; j++) {
    //                 row_buffer[j] = -row_buffer[j];
    //             }
    //         }
            
    //         // Copy sorted data back to local data array
    //         if (row_buffer) {
    //             int copy_size = (row_size < local_cols) ? row_size : local_cols;
    //             for (int j = 0; j < copy_size; j++) {
    //                 distrib->local_data[lr * local_cols + j] = row_buffer[j];
    //             }
    //             free(row_buffer);
    //         }
    //     }
        
    //     // Make sure all processes complete their row sort before proceeding to next row
    //     MPI_Barrier(distrib->row_comm);
    // }
    printf("Wating to be implemented for multiple processes per col.\n");
    for(int row_i = 0; row_i < local_rows; row_i++) {
        int global_row_idx = distrib->my_grid_row * local_rows + row_i;
        int is_descending = global_row_idx % 2; // Even rows asc, odd rows desc
        if(is_descending) {
            qsort(&(distrib->local_data[row_i * local_cols]), local_cols, sizeof(int), compare_descending);
        } else {
            qsort(&(distrib->local_data[row_i * local_cols]), local_cols, sizeof(int), compare_ascending);
        }

        
        
    }
    // Set a barrier to ensure all processes complete their row sort before proceeding to next row
    MPI_Barrier(distrib->row_comm);
    // Another loop to merge the sorted rows

    for(int row_i = 0; row_i < local_rows; row_i++) {
        int global_row_idx = distrib->my_grid_row * local_rows + row_i;
        int is_descending = global_row_idx % 2; // Even rows asc, odd rows desc
        if(is_descending) {
            // Merge the sorted rows in descending order
            global_sort(&(distrib->local_data[row_i * local_cols]), local_cols, distrib->row_comm, 1);
        } else {
            // Merge the sorted rows in ascending order
            global_sort(&(distrib->local_data[row_i * local_cols]), local_cols, distrib->row_comm, 0);
        }
    }



    return 0;
}