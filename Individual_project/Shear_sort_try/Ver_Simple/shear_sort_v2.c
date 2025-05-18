#include "shear_sort_v2.h"
#include "row_sort_v2.h"
#include "column_sort_v2.h"
#include "data_distribution_v2.h" // For gather_full_matrix
#include <stdio.h>
#include <stdlib.h>
#include <math.h>   // For ceil, log2
#include <string.h> // For strerror, potentially (not directly used here but good practice if errors expanded)
#include <mpi.h>

/**
 * Perform the Shear Sort algorithm (V2) on a distributed matrix.
 * Uses simplified row sort (v2) and column sort (v2).
 */
int perform_shear_sort_v2(DataDistribution *distrib, int pivot_strategy) {
    if (!distrib) {
        fprintf(stderr, "Error: NULL DataDistribution in perform_shear_sort_v2\n");
        return -1;
    }
    int N = distrib->N;
    int world_rank = distrib->world_rank;

    if (N == 0) return 0; // Nothing to sort for an empty matrix
    if (N < 0) {
        if (world_rank == 0) fprintf(stderr, "Error: Matrix size N cannot be negative (N=%d)\n", N);
        return -1;
    }

    // 1. Calculate number of phases: d = ceil(log2(N))
    // For N=1, log2(1)=0, d=0. Loop runs for phase_idx = 0. Row sort. Col sort skipped. Correct.
    // For N=2, log2(2)=1, d=1. Loop: phase 0 (row, col), phase 1 (row). Correct.
    int d = (N > 1) ? (int)ceil(log2((double)N)) : 0;

    // 2. Main iteration loop (Shear Sort requires log2(N) + 1 phases roughly)
    // The loop goes from phase_idx = 0 to d.
    // Total d+1 iterations for row sort.
    // Total d iterations for col sort.
    for (int phase_idx = 0; phase_idx <= d; phase_idx++) {
        if (distrib->world_rank == 0 && N > 32) { // Reduce verbosity for large N
            if (phase_idx % 5 == 0 || phase_idx == d) { // Print every 5 phases and the last one
                 printf("ShearSortV2: Phase %d/%d...\n", phase_idx, d);
            }
        }

        // 3. Row sort phase
        // perform_row_sort_v2 handles the snake-like direction based on global row index.
        int row_sort_result = perform_row_sort_v2(distrib, phase_idx);
        if (row_sort_result != 0) {
            fprintf(stderr, "Process %d: V2 Row sort failed in phase %d\n",
                    world_rank, phase_idx);
            return -1;
        }

        // 4. Column sort phase (always ascending)
        // Column sort is done in all phases except the very last one (phase d).
        if (phase_idx < d) {
            int col_sort_result = perform_column_sort_v2(distrib, phase_idx, pivot_strategy);
            if (col_sort_result != 0) {
                fprintf(stderr, "Process %d: V2 Column sort failed in phase %d after row sort\n",
                        world_rank, phase_idx);
                return -1;
            }
        }

        // Synchronize all processes before next phase if required by algorithm definition
        // (though individual sort routines might have their own internal syncs if they are complex collectives)
        MPI_Barrier(MPI_COMM_WORLD);
    }

    return 0;
}

/**
 * Verify if the matrix is sorted in snake-like order (V2).
 * Uses gather_full_matrix from data_distribution_v2.c.
 */
int verify_snake_order_v2(DataDistribution *distrib) {
    if (!distrib) {
        fprintf(stderr, "Error: NULL DataDistribution in verify_snake_order_v2\n");
        // Cannot determine rank to Bcast, so all might hang. Best to abort if possible.
        MPI_Abort(MPI_COMM_WORLD, 1); 
        return -1; // Should not be reached
    }
    int N = distrib->N;
    int world_rank = distrib->world_rank;
    int verification_status_on_root = 1; // Assume sorted, 0 if not, -1 if error during gather

    if (N == 0) { // An empty matrix is trivially sorted.
        if (world_rank == 0) printf("Verification: Empty matrix (N=0), considered sorted.\n");
        MPI_Bcast(&verification_status_on_root, 1, MPI_INT, 0, MPI_COMM_WORLD);
        return 1; 
    }
    if (N < 0) {
        if (world_rank == 0) fprintf(stderr, "Verification Error: Matrix size N cannot be negative (N=%d)\n", N);
        verification_status_on_root = 0; // Mark as not sorted / error
        MPI_Bcast(&verification_status_on_root, 1, MPI_INT, 0, MPI_COMM_WORLD);
        return 0;
    }

    int *global_matrix_on_root = NULL;

    if (world_rank == 0) {
        global_matrix_on_root = (int*)malloc(N * N * sizeof(int));
        if (global_matrix_on_root == NULL) {
            fprintf(stderr, "Root (verify_snake_order_v2): Memory allocation failed for verification buffer (N*N = %d*%d). Verification will fail.\n", N,N);
            verification_status_on_root = -1; // Error during setup
        }
    }

    // Broadcast root's allocation status (as part of verification_status_on_root)
    // If root failed malloc, other processes need to know to avoid MPI_ERR_BUFFER on gather.
    // A simpler way: root bcasts its `verification_status_on_root` now. If it's -1, all return.
    MPI_Bcast(&verification_status_on_root, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (verification_status_on_root == -1) { // Root malloc failed
        if (world_rank == 0 && global_matrix_on_root) free(global_matrix_on_root); // Should be NULL though
        return 0; // Indicate failure (not sorted because couldn't verify)
    }

    // All processes call gather_full_matrix (from data_distribution_v2.c)
    // Root passes the allocated buffer, others pass NULL (ignored by gather_full_matrix for non-roots).
    int gather_ret = gather_full_matrix(distrib, global_matrix_on_root);

    if (gather_ret != 0) {
        if (world_rank == 0) {
            fprintf(stderr, "Root (verify_snake_order_v2): gather_full_matrix failed (returned %d). Verification fails.\n", gather_ret);
        }
        verification_status_on_root = -1; // Gather error
    } else {
        // Gather succeeded. Only root performs the actual checks.
        if (world_rank == 0) {
            // Perform intra-row checks
            for (int i = 0; i < N; i++) { // Iterate global rows
                for (int j = 0; j < N - 1; j++) { // Iterate columns in the row
                    int current_val = global_matrix_on_root[i*N + j];
                    int next_val = global_matrix_on_root[i*N + j+1];
                    if (i % 2 == 0) { // Even global rows (0-indexed): should be ascending L-R
                        if (current_val > next_val) {
                            fprintf(stderr, "Verification Error: Row %d (even) not in ascending order at col %d vs %d: %d > %d\n",
                                    i, j, j+1, current_val, next_val);
                            verification_status_on_root = 0; break;
                        }
                    } else { // Odd global rows (0-indexed): should be descending L-R
                        if (current_val < next_val) {
                            fprintf(stderr, "Verification Error: Row %d (odd) not in descending order at col %d vs %d: %d < %d\n",
                                    i, j, j+1, current_val, next_val);
                            verification_status_on_root = 0; break;
                        }
                    }
                }
                if (!verification_status_on_root) break; // Exit outer loop if error found
            }

            // Perform inter-row (snake) connection checks if still considered sorted
            if (verification_status_on_root) {
                for (int i = 0; i < N - 1; i++) { // Iterate between row i and row i+1
                    int val_this_row_end, val_next_row_start;
                    if (i % 2 == 0) { // Current row i is even (L-R ascending), ends at N-1
                        val_this_row_end = global_matrix_on_root[i*N + (N-1)];
                        // Next row i+1 is odd (L-R descending), also "starts" check from its N-1 for snake logic based on shear sort spec
                        // No, for snake: end of even row <= start of odd row (which is element [i+1][N-1] if odd row sorted R-L, or [i+1][0] if L-R)
                        // The typical ShearSort snake means: R-end of even_row[i] <= R-end of odd_row[i+1] (if odd row R-L / descending L-R)
                        // OR L-end of odd_row[i] <= L-end of even_row[i+1]
                        // Standard snake: (even row i, L->R), (odd row i+1, R->L)
                        // global_matrix[i][N-1] <= global_matrix[i+1][N-1]
                        // If odd row (i+1) is L->R descending, then its largest is at [i+1][0]
                        // For our L->R ascending (even) and L->R descending (odd):
                        // End of even row (i*N + N-1) must be <= Start of odd row ((i+1)*N + 0) -- NO, this is for simple sorted matrix.
                        // For snake: last of even row <= last of odd row. OR first of odd row <= first of even row.
                        // Let's stick to the definition from original file for consistency of verification:
                        // Even row i (0,2,..) sorted L-R ascending. Connects to Odd row i+1 (1,3,..) sorted L-R descending.
                        // The connection is: element [i][N-1] <= element [i+1][N-1]. (Right ends)
                        val_this_row_end = global_matrix_on_root[i*N + (N-1)]; // Rightmost of current even row
                        val_next_row_start = global_matrix_on_root[(i+1)*N + (N-1)]; // Rightmost of next odd row
                    } else { // Current row i is odd (L-R descending), ends at 0
                        // Odd row i sorted L-R descending. Connects to Even row i+1 sorted L-R ascending.
                        // The connection is: element [i][0] <= element [i+1][0]. (Left ends)
                        val_this_row_end = global_matrix_on_root[i*N + 0];         // Leftmost of current odd row
                        val_next_row_start = global_matrix_on_root[(i+1)*N + 0];   // Leftmost of next even row
                    }

                    if (val_this_row_end > val_next_row_start) {
                         fprintf(stderr, "Verification Error: Snake-like order violated between row %d and row %d. Values: %d (end of row %d) > %d (start of row %d)\n",
                                 i, i+1, val_this_row_end, i, val_next_row_start, i+1);
                         verification_status_on_root = 0;
                         break;
                    }
                }
            }

            if (world_rank == 0) {
                if (verification_status_on_root == 1) {
                    printf("Verification V2 PASSED: Matrix is correctly sorted in snake-like order.\n");
                } else if (verification_status_on_root == 0) {
                    printf("Verification V2 FAILED: Matrix is NOT correctly sorted (see errors above).\n");
                } else { // -1
                    printf("Verification V2 FAILED: Error during matrix gather or setup.\n");
                }
            }
        } // end if (world_rank == 0) for verification logic
    } // end else (gather_ret == 0)

    if (world_rank == 0 && global_matrix_on_root != NULL) {
        free(global_matrix_on_root);
    }

    // All processes must have the final verification_status_on_root from root.
    MPI_Bcast(&verification_status_on_root, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Return 1 for sorted, 0 for not sorted/error for consistency with original verify. 
    // But verification_status_on_root now has -1 for gather error.
    // The main will check return value to print pass/fail.
    // Let's return 1 if sorted, 0 otherwise (which includes errors).
    return (verification_status_on_root == 1) ? 1 : 0;
} 