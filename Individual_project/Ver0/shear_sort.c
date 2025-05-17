#include "shear_sort.h"
#include "row_sort.h"
#include "column_sort.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/**
 * Perform the Shear Sort algorithm on a distributed matrix
 */
int perform_shear_sort(DataDistribution *distrib, int pivot_strategy) {
    int N = distrib->N;
    int world_rank = distrib->world_rank;
    
    // 1. Calculate number of phases: d = ceil(log2(N))
    int d = (int)ceil(log2((double)N));
    
    // 2. Main iteration loop (l from 1 to d+1 in the algorithm)
    for (int phase_idx = 0; phase_idx < d + 1; phase_idx++) {
        // 3. Row sort phase (both odd and even rows)
        // Note: algorithm uses 1-indexing, but we use 0-indexing
        // In the algorithm: odd rows (k=1,3,5...) are sorted ascending
        //                   even rows (k=2,4,6...) are sorted descending
        // In our code: even indices (0,2,4...) are sorted ascending
        //              odd indices (1,3,5...) are sorted descending
        // printf("Doing row sort in phase %d (pivot strategy %d) on process %d\n", 
        //       phase_idx, pivot_strategy, world_rank);
        int row_sort_result = perform_row_sort(distrib, phase_idx, pivot_strategy);
        if (row_sort_result != 0) {
            fprintf(stderr, "Process %d: Row sort failed in phase %d\n", 
                    world_rank, phase_idx);
            return -1;
        }
        // printf("Ending row sort in phase %d (pivot strategy %d) on process %d\n", phase_idx, pivot_strategy, world_rank);
        // 4. Column sort phase (only if not the last iteration)
        if (phase_idx < d) {
            int col_sort_result = perform_column_sort(distrib, phase_idx, pivot_strategy);
            if (col_sort_result != 0) {
                fprintf(stderr, "Process %d: Column sort failed in phase %d\n", 
                        world_rank, phase_idx);
                return -1;
            }
        }
        
        // Synchronize all processes before next phase
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    return 0;
}

/**
 * Gather the entire matrix on the root process for verification or output
 */
int gather_matrix(DataDistribution *distrib, int *output_global_matrix_on_root) {
    int N = distrib->N;
    int world_rank = distrib->world_rank;
    int world_size = distrib->world_size;
    int local_rows = distrib->local_rows;
    int local_cols = distrib->local_cols;
    
    // Check if root provided a valid buffer. This is crucial.
    if (world_rank == 0 && output_global_matrix_on_root == NULL) {
        fprintf(stderr, "Root (gather_matrix): output_global_matrix_on_root is NULL. This is a critical error. Aborting.\n");
        // All processes must be aware of this failure to prevent deadlock.
        MPI_Abort(MPI_COMM_WORLD, 1); // Abort all processes
        return -1; // Should not be reached if MPI_Abort is effective
    }
    
    // Create an MPI datatype for a matrix block
    MPI_Datatype block_type, resized_block_type;
    int sizes[2] = {N, N};              // Global matrix dimensions
    int subsizes[2] = {local_rows, local_cols};  // Local block dimensions
    int starts[2] = {0, 0};             // Starting offset
    
    // Create a datatype representing the local block layout in the global matrix
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &block_type);
    
    // Resize the datatype to represent the actual size of the local block
    MPI_Aint lb, extent;
    MPI_Type_get_extent(MPI_INT, &lb, &extent);
    MPI_Type_create_resized(block_type, 0, extent, &resized_block_type);
    MPI_Type_commit(&resized_block_type);
    
    // Calculate displacements and counts for gatherv
    int *displs = NULL;
    int *recvcounts = NULL;
    
    if (world_rank == 0) {
        displs = (int*)malloc(world_size * sizeof(int));
        recvcounts = (int*)malloc(world_size * sizeof(int));
        
        if (displs == NULL || recvcounts == NULL) {
            fprintf(stderr, "Root process: Memory allocation failed for gather parameters\n");
            free(displs);
            free(recvcounts);
            MPI_Type_free(&resized_block_type);
            return -1;
        }
        
        for (int i = 0; i < world_size; i++) {
            // Calculate grid coordinates for process i
            int proc_row = i / distrib->p_c;
            int proc_col = i % distrib->p_c;
            
            // Calculate displacement in terms of blocks
            displs[i] = proc_row * local_rows * N + proc_col * local_cols;
            
            // Each process sends one local block
            recvcounts[i] = 1;
        }
    }
    
    // Gather all local blocks to form the global matrix on root
    MPI_Gatherv(distrib->local_data, local_rows * local_cols, MPI_INT,
                output_global_matrix_on_root, // Root uses the provided buffer
                recvcounts, displs, resized_block_type,
                0, MPI_COMM_WORLD);
    
    // Clean up
    MPI_Type_free(&resized_block_type);
    if (world_rank == 0) {
        free(displs);
        free(recvcounts);
    }
    
    return 0;
}

/**
 * Verify if the matrix is sorted in snake-like order
 */
int verify_snake_order(DataDistribution *distrib) {
    int N = distrib->N;
    int world_rank = distrib->world_rank;
    int verification_result_on_root = 1; // Assume sorted initially

    int *global_matrix_for_verification = NULL;

    // Root allocates buffer. If fails, verification_result_on_root is set to 0.
    if (world_rank == 0) {
        global_matrix_for_verification = (int*)malloc(N * N * sizeof(int));
        if (global_matrix_for_verification == NULL) {
            fprintf(stderr, "Root (verify_snake_order): Memory allocation failed for verification buffer. Verification will fail.\n");
            verification_result_on_root = 0;
        }
    }

    // Broadcast root's allocation status. If root failed to allocate, all processes will skip gather and fail verification.
    // (More accurately, verification_result_on_root from root (0 if alloc failed, 1 if ok) is broadcast later).
    // The critical part is that gather_matrix must be called by all, or none if root can't provide buffer.
    // gather_matrix itself will MPI_Abort if root passes a NULL buffer.

    // So, if root failed malloc, global_matrix_for_verification is NULL.
    // The call to gather_matrix(distrib, NULL) from root will trigger MPI_Abort inside gather_matrix.
    // This is an acceptable way to handle root malloc failure for this collective operation chain.
    // If we wanted a gentler failure, root would Bcast a status, and all would skip gather_matrix.

    // ALL processes must call gather_matrix as it's collective.
    int gather_ret = gather_matrix(distrib, global_matrix_for_verification); // Non-root pass NULL, which is fine for their send part.

    if (gather_ret != 0) {
        // gather_matrix would have printed an error or aborted. 
        // If it returns non-zero (e.g. if root's buffer was NULL and it didn't abort but returned error),
        // then root should mark verification as failed.
        if (world_rank == 0) {
             fprintf(stderr, "Root (verify_snake_order): gather_matrix failed (returned %d). Verification fails.\n", gather_ret);
             verification_result_on_root = 0;
        }
        // If gather_matrix aborted, we might not reach here on all procs.
    } else {
        // Gather succeeded. Only root performs the actual checks if its buffer is valid.
        if (world_rank == 0) {
            if (global_matrix_for_verification == NULL) {
                // This case should ideally be caught by gather_matrix aborting if root passes NULL.
                // Or if root malloc failed AND gather_matrix somehow succeeded with a NULL buffer for root (very unlikely).
                fprintf(stderr, "Root (verify_snake_order): global_matrix_for_verification is NULL after supposedly successful gather. Verification fails.\n");
                verification_result_on_root = 0;
            } else {
                // Perform intra-row checks
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N - 1; j++) {
                        if (i % 2 == 0) { // Even rows (0-indexed): ascending L-R
                            if (global_matrix_for_verification[i*N + j] > global_matrix_for_verification[i*N + j+1]) {
                                fprintf(stderr, "Verification Error: Row %d (even) not in ascending order at col %d: %d > %d\n", 
                                        i, j, global_matrix_for_verification[i*N + j], global_matrix_for_verification[i*N + j+1]);
                                verification_result_on_root = 0; break;
                            }
                        } else { // Odd rows (0-indexed): descending L-R
                            if (global_matrix_for_verification[i*N + j] < global_matrix_for_verification[i*N + j+1]) {
                                fprintf(stderr, "Verification Error: Row %d (odd) not in descending order at col %d: %d < %d\n", 
                                        i, j, global_matrix_for_verification[i*N + j], global_matrix_for_verification[i*N + j+1]);
                                verification_result_on_root = 0; break;
                            }
                        }
                    }
                    if (!verification_result_on_root) break;
                }

                // Perform inter-row (snake) connection checks if still considered sorted
                if (verification_result_on_root) {
                    for (int i = 0; i < N - 1; i++) {
                        // Using the original logic for snake connections from the user's prior version of the code
                        int val_this_row, val_next_row;
                        if (i % 2 == 0) { // Current row i is even (L-R ascending)
                            val_this_row = global_matrix_for_verification[i*N + (N-1)];       // Rightmost element
                            // Next row i+1 is odd (L-R descending)
                            val_next_row = global_matrix_for_verification[(i+1)*N + (N-1)];   // Rightmost element (smallest for this L-R descending row)
                        } else { // Current row i is odd (L-R descending)
                            val_this_row = global_matrix_for_verification[i*N + 0];           // Leftmost element (largest for this L-R descending row)
                            // Next row i+1 is even (L-R ascending)
                            val_next_row = global_matrix_for_verification[(i+1)*N + 0];       // Leftmost element
                        }

                        if (val_this_row > val_next_row) {
                             fprintf(stderr, "Verification Error: Snake-like order violated between row %d and row %d. Values: %d > %d\n",
                                     i, i+1, val_this_row, val_next_row);
                             verification_result_on_root = 0;
                             break;
                        }
                    }
                }

                if (verification_result_on_root) {
                    printf("Verification passed: Matrix is correctly sorted in snake-like order\n");
                } else {
                    printf("Verification failed: Matrix is NOT correctly sorted (see errors above or gather failed).\n");
                }
            } // end else (global_matrix_for_verification != NULL)
        } // end if (world_rank == 0) for verification logic
    } // end else (gather_ret == 0)

    if (world_rank == 0 && global_matrix_for_verification != NULL) {
        free(global_matrix_for_verification);
    }

    // All processes must participate in this Bcast.
    MPI_Bcast(&verification_result_on_root, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    return verification_result_on_root; 
}