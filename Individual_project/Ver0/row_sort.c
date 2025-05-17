#include "row_sort.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/**
 * Merge two arrays in ascending order
 */
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

/**
 * Merge two arrays in descending order
 */
void merge_descending(int *v1, int n1, int *v2, int n2, int *result) {
    int i = 0, j = 0, k = 0;
    while (i < n1 && j < n2) {
        if (v1[i] >= v2[j]) {
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

/**
 * Helper function to get the median of an array
 */
int get_median(int *arr, int n) {
    // Simple insertion sort for small arrays
    int *temp = malloc(n * sizeof(int));
    memcpy(temp, arr, n * sizeof(int));
    
    for (int i = 1; i < n; i++) {
        int key = temp[i];
        int j = i - 1;
        while (j >= 0 && temp[j] > key) {
            temp[j + 1] = temp[j];
            j = j - 1;
        }
        temp[j + 1] = key;
    }
    
    int median = temp[n/2];
    free(temp);
    return median;
}

/**
 * Find the index of the first element larger than the pivot value
 * Used for partitioning
 */
int get_larger_index(int *elements, int n, int pivot_val) {
    if (n <= 0 || elements == NULL) return 0;
    
    int i;
    for (i = 0; i < n; i++) {
        if (elements[i] > pivot_val) {
            break;
        }
    }
    return i;
}

/**
 * Select pivot using mean of medians strategy
 */
int select_pivot_mean_median(int *elements, int n, MPI_Comm communicator) {
    int rank, size;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    int local_median = 0; 
    if (n > 0 && elements != NULL) {
        local_median = get_median(elements, n); 
    }

    int *all_medians = NULL;
    if (rank == ROOT && size > 0) { 
        all_medians = (int *)malloc(size * sizeof(int));
    }

    MPI_Gather(&local_median, 1, MPI_INT, all_medians, 1, MPI_INT, ROOT, communicator);

    int pivot_val = 0; 
    if (rank == ROOT) {
        long long sum_of_medians = 0;
        for (int i = 0; i < size; ++i) {
            sum_of_medians += all_medians[i];
        }
        pivot_val = (int)(sum_of_medians / size);
        free(all_medians);
    }

    MPI_Bcast(&pivot_val, 1, MPI_INT, ROOT, communicator);
    return get_larger_index(elements, n, pivot_val);
}

/**
 * Select pivot based on strategy
 */
int select_pivot(int pivot_strategy, int *elements, int n, MPI_Comm communicator) {
    // This function is now only expected to return THE PIVOT VALUE.
    // The select_pivot_mean_median needs to be refactored or this function adapted.
    // For now, to ensure select_pivot_and_partition gets a value, let's replicate pivot value calculation here.
    int rank_in_comm, size_in_comm;
    MPI_Comm_rank(communicator, &rank_in_comm);
    MPI_Comm_size(communicator, &size_in_comm);

    int local_median_for_pivot = 0; 
    if (n > 0 && elements != NULL) {
        local_median_for_pivot = get_median(elements, n); 
    }

    int *all_medians_for_pivot = NULL;
    if (rank_in_comm == ROOT && size_in_comm > 0) { 
        all_medians_for_pivot = (int *)malloc(size_in_comm * sizeof(int));
        if (!all_medians_for_pivot && size_in_comm > 0) { MPI_Abort(communicator, 1); return 0;}
    }

    MPI_Gather(&local_median_for_pivot, 1, MPI_INT, all_medians_for_pivot, 1, MPI_INT, ROOT, communicator);

    int determined_pivot_val = 0; 
    if (rank_in_comm == ROOT) {
        long long sum_of_medians_val = 0;
        if (size_in_comm > 0) {
            for (int i = 0; i < size_in_comm; ++i) {
                sum_of_medians_val += all_medians_for_pivot[i];
            }
            determined_pivot_val = (int)(sum_of_medians_val / size_in_comm);
        } else {
            determined_pivot_val = 0; 
        }
        if (all_medians_for_pivot) free(all_medians_for_pivot);
    }

    MPI_Bcast(&determined_pivot_val, 1, MPI_INT, ROOT, communicator);
    return determined_pivot_val;
}

/**
 * Select pivot and partition locally
 * Modified to support both ascending and descending sorting
 */
int select_pivot_and_partition(int **elements_ptr, int n, MPI_Comm communicator,
                               int pivot_strategy, int *num_small_local, int *num_large_local,
                               int sort_direction) {
    if (n == 0) {
        *num_small_local = 0;
        *num_large_local = 0;
        return 0; // No pivot index needed, or count is 0
    }

    // 1. Determine the actual pivot value using the (now modified) select_pivot function.
    int pivot_val = select_pivot(pivot_strategy, *elements_ptr, n, communicator);

    // 2. Partition the local *elements_ptr array based on pivot_val and sort_direction.
    int *temp_copy = (int*)malloc(n * sizeof(int));
    if (!temp_copy) { 
        fprintf(stderr, "Process (select_pivot_and_partition): Malloc failed for temp_copy\n");
        MPI_Abort(communicator, 1); 
        return -1; // Should not be reached
    }
    memcpy(temp_copy, *elements_ptr, n * sizeof(int));

    int current_partition_idx = 0; // This will be the count of elements in the "first" group after partitioning
    int k = 0; // writer index for *elements_ptr

    if (sort_direction == SORT_ASCENDING) {
        // "Small" group (elements <= pivot_val) comes first, then "Large" group (elements > pivot_val)
        for (int j = 0; j < n; j++) { // First pass for elements <= pivot_val
            if (temp_copy[j] <= pivot_val) {
                (*elements_ptr)[k++] = temp_copy[j];
            }
        }
        current_partition_idx = k; // Count of elements in the "small" group
        for (int j = 0; j < n; j++) { // Second pass for elements > pivot_val
            if (temp_copy[j] > pivot_val) {
                (*elements_ptr)[k++] = temp_copy[j];
            }
        }
    } else { // SORT_DESCENDING
        // "Small" group (elements >= pivot_val) comes first, then "Large" group (elements < pivot_val)
        for (int j = 0; j < n; j++) { // First pass for elements >= pivot_val
            if (temp_copy[j] >= pivot_val) {
                (*elements_ptr)[k++] = temp_copy[j];
            }
        }
        current_partition_idx = k; // Count of elements in the "small" (i.e., >= pivot) group
        for (int j = 0; j < n; j++) { // Second pass for elements < pivot_val
            if (temp_copy[j] < pivot_val) {
                (*elements_ptr)[k++] = temp_copy[j];
            }
        }
    }
    
    free(temp_copy);

    *num_small_local = current_partition_idx;
    *num_large_local = n - current_partition_idx;
    
    // The return value of this function isn't strictly used by global_sort in its current form,
    // as num_small_local and num_large_local are output parameters.
    return current_partition_idx; 
}

/**
 * Determine process groups, partners, and send/keep buffers
 * Modified to support both ascending and descending sorting
 */
void determine_process_groups(int my_rank, int num_procs, int num_small_local, int num_large_local,
                              int **elements_ptr, int *group, int *partner_rank,
                              int **send_buf_ptr, int *send_count,
                              int **keep_buf_start_ptr, int *keep_count,
                              int sort_direction) {
    *group = (my_rank < num_procs / 2) ? 0 : 1;

    // For ascending sort, low group keeps smaller elements, high group keeps larger elements
    // For descending sort, low group keeps larger elements, high group keeps smaller elements
    
    if (*group == 0) {  // Low group
        *partner_rank = my_rank + num_procs / 2;
        
        if (sort_direction == SORT_ASCENDING) {
            // Low group keeps small elements, sends large elements
            *send_buf_ptr = (*elements_ptr) + num_small_local;
            *send_count = num_large_local;
            *keep_buf_start_ptr = *elements_ptr;
            *keep_count = num_small_local;
        } else {
            // For descending, low group keeps large elements, sends small elements
            *send_buf_ptr = (*elements_ptr) + num_small_local;
            *send_count = num_large_local;
            *keep_buf_start_ptr = *elements_ptr;
            *keep_count = num_small_local;
        }
    } else {  // High group
        *partner_rank = my_rank - num_procs / 2;
        
        if (sort_direction == SORT_ASCENDING) {
            // High group keeps large elements, sends small elements
            *send_buf_ptr = *elements_ptr;
            *send_count = num_small_local;
            *keep_buf_start_ptr = (*elements_ptr) + num_small_local;
            *keep_count = num_large_local;
        } else {
            // For descending, high group keeps small elements, sends large elements
            *send_buf_ptr = *elements_ptr;
            *send_count = num_small_local;
            *keep_buf_start_ptr = (*elements_ptr) + num_small_local;
            *keep_count = num_large_local;
        }
    }
}

/**
 * Exchange data with partner
 */
int* exchange_data(int *send_buf_ptr, int send_count, int partner_rank,
                   MPI_Comm communicator, int my_rank, int *recv_count_ptr) {
    int recv_count;
    // Exchange send_count with partner to determine recv_count
    MPI_Sendrecv(&send_count, 1, MPI_INT, partner_rank, 0,
                 &recv_count, 1, MPI_INT, partner_rank, 0,
                 communicator, MPI_STATUS_IGNORE);

    *recv_count_ptr = recv_count;
    int *recv_buf = NULL;
    
    // Allocate buffer for received data if recv_count > 0
    if (recv_count > 0) {
        recv_buf = malloc(recv_count * sizeof(int));
        if (!recv_buf) {
            fprintf(stderr, "Memory allocation failed in exchange_data\n");
            MPI_Abort(communicator, 1);
        }
    }
    
    // Exchange actual data
    MPI_Sendrecv(send_buf_ptr, send_count, MPI_INT, partner_rank, 1,
                 recv_buf, recv_count, MPI_INT, partner_rank, 1,
                 communicator, MPI_STATUS_IGNORE);

    return recv_buf;
}

/**
 * Merge local kept data with received data
 * Modified to support both ascending and descending sorting
 */
int* merge_data(int *keep_buf_start_ptr, int keep_count, int *recv_buf, int recv_count,
                int my_rank, MPI_Comm communicator, int *new_n_ptr, int sort_direction) {
    int new_n = keep_count + recv_count;
    *new_n_ptr = new_n;
    int *new_elements_array = NULL;

    if (new_n > 0) {
        // Allocate new array for merged data
        new_elements_array = malloc(new_n * sizeof(int));
        if (!new_elements_array) {
            fprintf(stderr, "Memory allocation error in merge_data\n");
            MPI_Abort(communicator, 1);
            return NULL;
        }
        
        if (keep_count > 0 && recv_count > 0) {
            // Merge two sorted arrays
            int *temp_array = malloc(keep_count * sizeof(int));
            if (!temp_array) {
                fprintf(stderr, "Memory allocation error in merge_data\n");
                MPI_Abort(communicator, 1);
                free(new_elements_array);
                return NULL;
            }
            
            memcpy(temp_array, keep_buf_start_ptr, keep_count * sizeof(int));
            
            // Use appropriate merge function based on sort direction
            if (sort_direction == SORT_ASCENDING) {
                merge_ascending(temp_array, keep_count, recv_buf, recv_count, new_elements_array);
            } else {
                merge_descending(temp_array, keep_count, recv_buf, recv_count, new_elements_array);
            }
            
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
 * Split communicator and recurse
 */
int recursive_sort(int **elements_ptr, int n, MPI_Comm communicator,
                   int group, int my_rank, int pivot_strategy, int sort_direction) {
    // Create new communicator for the group
    MPI_Comm new_comm;
    MPI_Comm_split(communicator, group, my_rank, &new_comm);

    // Sort in the new communicator
    n = global_sort(elements_ptr, n, new_comm, pivot_strategy, sort_direction);
    
    // Free the new communicator
    MPI_Comm_free(&new_comm);
    return n;
}

/**
 * Main parallel quicksort
 * Modified to support both ascending and descending sorting
 */
int global_sort(int **elements_ptr, int n, MPI_Comm communicator, 
               int pivot_strategy, int sort_direction) {
    int my_rank, num_procs;
    MPI_Comm_rank(communicator, &my_rank);
    MPI_Comm_size(communicator, &num_procs);

    // Base case: if only one process in the communicator, sort locally and return.
    // If n <= 1 on this single process, the local sort will do nothing or handle it correctly.
    if (num_procs <= 1) { 
        if (n > 1) { // Only sort if there's more than one element locally
            // Simple insertion sort for small arrays
            for (int i = 1; i < n; i++) {
                int key = (*elements_ptr)[i];
                int j = i - 1;
                
                if (sort_direction == SORT_ASCENDING) {
                    while (j >= 0 && (*elements_ptr)[j] > key) {
                        (*elements_ptr)[j + 1] = (*elements_ptr)[j];
                        j--;
                    }
                } else { // SORT_DESCENDING
                    while (j >= 0 && (*elements_ptr)[j] < key) {
                        (*elements_ptr)[j + 1] = (*elements_ptr)[j];
                        j--;
                    }
                }
                (*elements_ptr)[j + 1] = key;
            }
        }
        return n; // Return current number of elements
    }

    // If n is 0 for all processes in this communicator, we might also be able to return early.
    // However, the recursive structure should handle this: if n=0, num_small_local and num_large_local will be 0,
    // send_count will be 0, recv_count will be 0, new_n will be 0, and recursion continues with n=0.
    // This seems acceptable. The primary base case is num_procs <= 1.

    // 1. Select pivot and partition locally
    int num_small_local, num_large_local;
    select_pivot_and_partition(elements_ptr, n, communicator,
                              pivot_strategy, &num_small_local, &num_large_local,
                              sort_direction);

    // 2. Determine groups and buffers
    int group, partner_rank;
    int *send_buf_ptr, *keep_buf_start_ptr;
    int send_count, keep_count;
    determine_process_groups(my_rank, num_procs, num_small_local, num_large_local,
                           elements_ptr, &group, &partner_rank,
                           &send_buf_ptr, &send_count,
                           &keep_buf_start_ptr, &keep_count,
                           sort_direction);

    // 3. Exchange data
    int recv_count;
    int *recv_buf = exchange_data(send_buf_ptr, send_count, partner_rank,
                                communicator, my_rank, &recv_count);

    // 4. Merge data
    int new_n;
    int *new_elements_array = merge_data(keep_buf_start_ptr, keep_count,
                                      recv_buf, recv_count,
                                      my_rank, communicator, &new_n,
                                      sort_direction);

    // Free old elements array and assign new one
    free(*elements_ptr);
    *elements_ptr = new_elements_array;
    n = new_n;

    // Free received buffer if it was allocated
    if (recv_buf) {
        free(recv_buf);
    }

    // 5. Recurse on new communicator
    n = recursive_sort(elements_ptr, n, communicator, group, my_rank, 
                    pivot_strategy, sort_direction);

    return n;
}

/**
 * Perform the row sorting phase of shearsort
 */
int perform_row_sort(DataDistribution *distrib, int phase_idx, int pivot_strategy) {
    int local_rows = distrib->local_rows;
    int local_cols = distrib->local_cols;
    int N = distrib->N;
    
    // For each local row
    for (int local_row = 0; local_row < local_rows; local_row++) {
        // Calculate global row number
        int global_row = distrib->my_grid_row * local_rows + local_row;
        
        // Determine sort direction based on row parity
        // Even rows (0-indexed) are sorted in ascending order
        // Odd rows (0-indexed) are sorted in descending order
        int sort_direction = (global_row % 2 == 0) ? SORT_ASCENDING : SORT_DESCENDING;
        
        // Extract the row data
        int *row_segment = (int*)malloc(local_cols * sizeof(int));
        if (!row_segment) {
            fprintf(stderr, "Process %d: Memory allocation failed for row segment\n", 
                    distrib->world_rank);
            return -1;
        }
        
        // Copy row data to the segment buffer
        memcpy(row_segment, &(distrib->local_data[local_row * local_cols]), 
               local_cols * sizeof(int));
        
        // Sort the row using global_sort
        int num_elements = local_cols;
        int *sorted_segment = row_segment;  // Start with the extracted segment
        
        // Perform global sort on this row using row communicator
        // global_sort will free the memory pointed to by sorted_segment (i.e., row_segment's content)
        // and allocate new memory for the result, updating sorted_segment to point to it.
        // num_elements will be the count of items this process holds for the sorted row.
        num_elements = global_sort(&sorted_segment, num_elements, 
                                 distrib->row_comm, pivot_strategy, sort_direction);
        
        // At this point, sorted_segment contains 'num_elements' sorted items for this process's
        // part of the current global row. The sum of 'num_elements' across all procs in row_comm is N.
        // We need to redistribute these so each process gets 'local_cols' elements back, forming
        // its correct part of the globally sorted row.

        int my_rank_in_row_comm, num_procs_in_row_comm;
        MPI_Comm_rank(distrib->row_comm, &my_rank_in_row_comm);
        MPI_Comm_size(distrib->row_comm, &num_procs_in_row_comm);

        // 1. Get all counts of elements each process in row_comm holds after global_sort
        int *all_proc_element_counts = (int*)malloc(num_procs_in_row_comm * sizeof(int));
        if (!all_proc_element_counts) {
            fprintf(stderr, "Process %d: Memory allocation failed for all_proc_element_counts\n", distrib->world_rank);
            if (sorted_segment) free(sorted_segment); // sorted_segment was allocated by global_sort
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }

        MPI_Allgather(&num_elements, 1, MPI_INT, 
                      all_proc_element_counts, 1, MPI_INT, distrib->row_comm);

        // 2. Calculate displacements for MPI_Allgatherv to reconstruct the full sorted row
        int *displacements = (int*)malloc(num_procs_in_row_comm * sizeof(int));
        if (!displacements) {
            fprintf(stderr, "Process %d: Memory allocation failed for displacements\n", distrib->world_rank);
            free(all_proc_element_counts);
            if (sorted_segment) free(sorted_segment);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }

        displacements[0] = 0;
        for (int i = 1; i < num_procs_in_row_comm; i++) {
            displacements[i] = displacements[i-1] + all_proc_element_counts[i-1];
        }
        
        // Optional: Verify total elements. Sum of all_proc_element_counts should be N.
        // int total_row_elements_check = 0;
        // for(int i=0; i<num_procs_in_row_comm; ++i) total_row_elements_check += all_proc_element_counts[i];
        // if (total_row_elements_check != N) {
        //    fprintf(stderr, "Process %d (Global Row %d): MISMATCH - total elements after sort in row_comm (%d) != N (%d)\n", 
        //            distrib->world_rank, global_row, total_row_elements_check, N);
        // }

        // 3. Gather all pieces of the sorted row from all processes in row_comm into a full_sorted_row buffer.
        // The size of the full_sorted_row is N (distrib->N).
        int *full_sorted_row = (int*)malloc(N * sizeof(int));
        if (!full_sorted_row) {
            fprintf(stderr, "Process %d: Memory allocation failed for full_sorted_row\n", distrib->world_rank);
            free(all_proc_element_counts);
            free(displacements);
            if (sorted_segment) free(sorted_segment);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }

        MPI_Allgatherv(sorted_segment, num_elements, MPI_INT,             // Current process sends its piece
                       full_sorted_row, all_proc_element_counts, displacements, MPI_INT, // Receives all pieces
                       distrib->row_comm);

        // 4. Copy the correct segment (local_cols elements) from full_sorted_row back to distrib->local_data.
        // The starting offset for this process in full_sorted_row is its rank in row_comm * local_cols.
        int start_offset_in_full_row = my_rank_in_row_comm * local_cols;

        if (start_offset_in_full_row + local_cols <= N) {
            memcpy(&(distrib->local_data[local_row * local_cols]), 
                   &full_sorted_row[start_offset_in_full_row], 
                   local_cols * sizeof(int));
        } else {
            // This should not happen if N is divisible by num_procs_in_row_comm (p_c)
            // and local_cols is correctly N / p_c.
            fprintf(stderr, "Process %d (Rank in RowComm %d): Error - Data copy bounds error. Offset %d + local_cols %d > N %d for global_row %d\n",
                    distrib->world_rank, my_rank_in_row_comm, start_offset_in_full_row, local_cols, N, global_row);
            free(full_sorted_row);
            free(all_proc_element_counts);
            free(displacements);
            if (sorted_segment) free(sorted_segment);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }
        
        // Free the buffers used in this iteration
        if (sorted_segment != NULL) { // sorted_segment was allocated by global_sort
            free(sorted_segment);
            sorted_segment = NULL; // Avoid double free if loop continues/errors
        }
        free(all_proc_element_counts);
        free(displacements);
        free(full_sorted_row);

        // The old error check is removed as we now handle variable num_elements by redistributing.
        // if (num_elements == local_cols) { ...
    }
    
    // Synchronize after all rows are sorted by all processes in MPI_COMM_WORLD
    MPI_Barrier(MPI_COMM_WORLD);
    
    return 0;
}