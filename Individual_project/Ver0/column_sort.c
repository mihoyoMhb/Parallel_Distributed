#include "column_sort.h"
#include "row_sort.h"  // For global_sort and sort direction constants
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Perform the column sorting phase of shearsort
 * Column sorting is always in ascending order
 */
int perform_column_sort(DataDistribution *distrib, int phase_idx, int pivot_strategy) {
    int local_rows = distrib->local_rows;
    int local_cols = distrib->local_cols;
    int N = distrib->N;
    
    // For each local column
    for (int local_col = 0; local_col < local_cols; local_col++) {
        // Calculate global column number
        int global_col = distrib->my_grid_col * local_cols + local_col;
        
        // Extract the column data (non-contiguous in memory)
        int *col_segment = (int*)malloc(local_rows * sizeof(int));
        if (!col_segment) {
            fprintf(stderr, "Process %d: Memory allocation failed for column segment\n", 
                    distrib->world_rank);
            return -1;
        }
        
        // Copy column data to the segment buffer (non-contiguous elements with stride local_cols)
        for (int i = 0; i < local_rows; i++) {
            col_segment[i] = distrib->local_data[i * local_cols + local_col];
        }
        
        // Sort the column using global_sort (always ascending for columns)
        int num_elements = local_rows;
        int *sorted_segment = col_segment;  // Start with the extracted segment
        
        // Perform global sort on this column using column communicator
        // global_sort will free the memory pointed to by sorted_segment (i.e., col_segment's content)
        // and allocate new memory for the result, updating sorted_segment to point to it.
        // num_elements will be the count of items this process holds for the sorted column.
        num_elements = global_sort(&sorted_segment, num_elements, 
                                  distrib->col_comm, pivot_strategy, SORT_ASCENDING);
        
        // At this point, sorted_segment contains 'num_elements' sorted items for this process's
        // part of the current global column. The sum of 'num_elements' across all procs in col_comm is N.
        // We need to redistribute these so each process gets 'local_rows' elements back, forming
        // its correct part of the globally sorted column.

        int my_rank_in_col_comm, num_procs_in_col_comm;
        MPI_Comm_rank(distrib->col_comm, &my_rank_in_col_comm);
        MPI_Comm_size(distrib->col_comm, &num_procs_in_col_comm);

        // 1. Get all counts of elements each process in col_comm holds after global_sort
        int *all_proc_element_counts = (int*)malloc(num_procs_in_col_comm * sizeof(int));
        if (!all_proc_element_counts) {
            fprintf(stderr, "Process %d: Memory allocation failed for all_proc_element_counts (column sort)\n", distrib->world_rank);
            if (sorted_segment) free(sorted_segment);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }

        MPI_Allgather(&num_elements, 1, MPI_INT, 
                      all_proc_element_counts, 1, MPI_INT, distrib->col_comm);

        // 2. Calculate displacements for MPI_Allgatherv to reconstruct the full sorted column
        int *displacements = (int*)malloc(num_procs_in_col_comm * sizeof(int));
        if (!displacements) {
            fprintf(stderr, "Process %d: Memory allocation failed for displacements (column sort)\n", distrib->world_rank);
            free(all_proc_element_counts);
            if (sorted_segment) free(sorted_segment);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }

        displacements[0] = 0;
        for (int i = 1; i < num_procs_in_col_comm; i++) {
            displacements[i] = displacements[i-1] + all_proc_element_counts[i-1];
        }

        // Optional: Verify total elements. Sum of all_proc_element_counts should be N.
        // int total_col_elements_check = 0;
        // for(int i=0; i<num_procs_in_col_comm; ++i) total_col_elements_check += all_proc_element_counts[i];
        // if (total_col_elements_check != N) {
        //    fprintf(stderr, "Process %d (Global Col %d): MISMATCH - total elements after sort in col_comm (%d) != N (%d)\n", 
        //            distrib->world_rank, global_col, total_col_elements_check, N);
        // }

        // 3. Gather all pieces of the sorted column from all processes in col_comm into a full_sorted_col buffer.
        // The size of the full_sorted_col is N (distrib->N, which is the global number of rows).
        int *full_sorted_col = (int*)malloc(N * sizeof(int));
        if (!full_sorted_col) {
            fprintf(stderr, "Process %d: Memory allocation failed for full_sorted_col (column sort)\n", distrib->world_rank);
            free(all_proc_element_counts);
            free(displacements);
            if (sorted_segment) free(sorted_segment);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }

        MPI_Allgatherv(sorted_segment, num_elements, MPI_INT,             // Current process sends its piece
                       full_sorted_col, all_proc_element_counts, displacements, MPI_INT, // Receives all pieces
                       distrib->col_comm);

        // 4. Copy the correct segment (local_rows elements) from full_sorted_col back to distrib->local_data.
        // The starting offset for this process in full_sorted_col is its rank in col_comm * local_rows.
        int start_offset_in_full_col = my_rank_in_col_comm * local_rows;

        if (start_offset_in_full_col + local_rows <= N) {
            for (int i = 0; i < local_rows; i++) {
                distrib->local_data[i * local_cols + local_col] = full_sorted_col[start_offset_in_full_col + i];
            }
        } else {
            fprintf(stderr, "Process %d (Rank in ColComm %d): Error - Column data copy bounds. Offset %d + local_rows %d > N %d for global_col %d\n",
                    distrib->world_rank, my_rank_in_col_comm, start_offset_in_full_col, local_rows, N, global_col);
            free(full_sorted_col);
            free(all_proc_element_counts);
            free(displacements);
            if (sorted_segment) free(sorted_segment);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }
        
        // Free the buffers used in this iteration
        if (sorted_segment != NULL) { // sorted_segment was allocated by global_sort
            free(sorted_segment);
            sorted_segment = NULL; 
        }
        free(all_proc_element_counts);
        free(displacements);
        free(full_sorted_col);
    }
    
    // Synchronize after all columns are sorted
    MPI_Barrier(MPI_COMM_WORLD);
    
    return 0;
}