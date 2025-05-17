#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <string.h> // For strerror
#include <errno.h>  // For errno
#include "data_distribution.h"
#include "row_sort.h"
#include "column_sort.h"
#include "shear_sort.h"

#define DEFAULT_MATRIX_SIZE 16 // This will be overridden by CSV if present
#define CSV_FILENAME "matrix_data.csv"

// Print a matrix (for debugging purposes)
void print_global_matrix_on_root(int *matrix, int N, int rank) {
    if (rank == 0) {
        printf("Global Matrix (%dx%d) on Root:\n", N, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%6d ", matrix[i * N + j]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
    int N = DEFAULT_MATRIX_SIZE;
    int rank, size;
    double start_time, end_time;
    const char* csv_filename = CSV_FILENAME;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Optional: Allow overriding CSV filename from command line
    // if (argc > 1) { csv_filename = argv[1]; }

    if (rank == 0) {
        FILE *fp = fopen(csv_filename, "r");
        if (fp == NULL) {
            fprintf(stderr, "Root: Error opening CSV file %s: %s. Using default N=%d\n", 
                    csv_filename, strerror(errno), DEFAULT_MATRIX_SIZE);
            // N remains DEFAULT_MATRIX_SIZE if file can't be opened
        } else {
            int csv_rows, csv_cols;
            if (fscanf(fp, "%d,%d", &csv_rows, &csv_cols) == 2) {
                if (csv_rows != csv_cols) {
                    fprintf(stderr, "Root: CSV matrix dimensions %dx%d are not square. Using first dimension %d.\n", 
                            csv_rows, csv_cols, csv_rows);
                    // Or handle error more strictly
                }
                N = csv_rows; // Use the row count as N for a square matrix
                printf("Root: Read N=%d from %s\n", N, csv_filename);
            } else {
                fprintf(stderr, "Root: Error reading dimensions from %s. Using default N=%d\n", 
                        csv_filename, DEFAULT_MATRIX_SIZE);
            }
            fclose(fp);
        }
    }

    // Broadcast N from root to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (N <= 0) {
        if (rank == 0) fprintf(stderr, "Error: Matrix size N must be positive. N=%d\n", N);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Initialize data distribution
    DataDistribution distrib;
    if (init_distribution(N, &distrib) != 0) {
        fprintf(stderr, "Process %d: Failed to initialize distribution for N=%d\n", rank, N);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (setup_processor_grid(&distrib) != 0) {
        fprintf(stderr, "Process %d: Failed to setup processor grid\n", rank);
        // Check for p_r or p_c being zero if N is too small for the number of processes
        if (distrib.p_r == 0 || distrib.p_c == 0) {
             if (rank==0) fprintf(stderr, "Error: Too many processes for matrix size N=%d. p_r=%d, p_c=%d\n", N, distrib.p_r, distrib.p_c);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (distrib.local_rows == 0 || distrib.local_cols == 0) {
        if (rank == 0) {
            fprintf(stderr, "Error: Calculated local block dimensions are zero (local_rows=%d, local_cols=%d). \n", 
                    distrib.local_rows, distrib.local_cols);
            fprintf(stderr, "This can happen if N (%d) is too small for the processor grid (%dx%d).\n", 
                    N, distrib.p_r, distrib.p_c);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (allocate_local_block(&distrib) != 0) {
        fprintf(stderr, "Process %d: Failed to allocate local block\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // Read data from CSV (root reads) and distribute to all processes
    if (read_csv_and_distribute_data(csv_filename, &distrib) != 0) {
        fprintf(stderr, "Process %d: Failed to read/distribute data from CSV %s\n", rank, csv_filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Print initial local block for each process
    print_local_block(&distrib, "Initial local block");
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize before printing global (if done)

    int *initial_full_matrix = NULL;
    // This whole block for gathering and printing initial global matrix is conditional on N <= 32
    if (N <= 32) {
        if (rank == 0) {
            // Root allocates the buffer for gathering
            initial_full_matrix = (int*)malloc(N * N * sizeof(int));
            if (initial_full_matrix == NULL) {
                fprintf(stderr, "Root: Failed to allocate memory for initial_full_matrix. Printing of gathered global matrix will be skipped.\n");
                // initial_full_matrix remains NULL. gather_matrix in shear_sort.c might handle this or fail.
                // For safety, if malloc fails, root should ideally not proceed with a gather that needs this buffer,
                // or ensure gather_matrix is robust to a NULL buffer on root (e.g. by not attempting the MPI_Gatherv).
                // However, to fix the immediate hang, we ensure all call gather_matrix.
                // The robustness of gather_matrix to a NULL buffer on root is a separate concern (Issue #2).
            }
        }

        // ALL processes must call gather_matrix because it uses MPI_Gatherv (a collective operation).
        // The 'initial_full_matrix' argument is the receive buffer for the root process.
        // For non-root processes, this argument is effectively ignored by MPI_Gatherv's recvbuf parameter.
        // So, non-root processes passing their local NULL initial_full_matrix is fine.
        int gather_res = gather_matrix(&distrib, initial_full_matrix);
        if (gather_res == 0) {
            if (rank == 0 && initial_full_matrix != NULL) {
                // Only root prints, and only if its buffer is valid
                // print_global_matrix_on_root(initial_full_matrix, N, rank);
            }
        } else {
            if (rank == 0) {
                fprintf(stderr, "Root: gather_matrix failed for initial data (returned %d).\n", gather_res);
            }
            // Other ranks could also log error if gather_matrix returns non-zero for them.
        }
    } // End of if (N <= 32)
    
    // User debug print - ensure ordered printing for clarity
    // for (int i = 0; i < size; ++i) {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (rank == i) {
    //         printf("Process %d: Print for debugging after initial gather block (N=%d)\n", rank, N);
    //     }
    // }
    // MPI_Barrier(MPI_COMM_WORLD); // Ensure all prints are done before proceeding

    // MPI_Barrier(MPI_COMM_WORLD); // This was the original barrier where the hang was reported.
                                 // It should now be fine after the gather_matrix fix.
    if (rank == 0) {
        printf("Starting Shear Sort...\n");
    }
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes before timing and starting sort

    start_time = MPI_Wtime();
    int pivot_strategy = 0; 
    int result = perform_shear_sort(&distrib, pivot_strategy);
    
    end_time = MPI_Wtime();
    
    if (result != 0) {
        fprintf(stderr, "Process %d: Shear sort failed\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (rank == 0) {
        printf("Shear sort on %dx%d matrix completed in %f seconds\n", N, N, end_time - start_time);
    }

    // Print final local block for each process
    print_local_block(&distrib, "Final local block (after sort)");
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize before printing global (if done)

    int *final_full_matrix = NULL;
    // This whole block for gathering and printing final global matrix is conditional on N <= 32
    if (N <= 32) {
        if (rank == 0) {
            final_full_matrix = (int*)malloc(N * N * sizeof(int));
            if (final_full_matrix == NULL) {
                fprintf(stderr, "Root: Failed to allocate memory for final_full_matrix. Printing of gathered global matrix will be skipped.\n");
            }
        }

        // ALL processes must call gather_matrix for the final matrix as well.
        int gather_res_final = gather_matrix(&distrib, final_full_matrix);
        if (gather_res_final == 0) {
            if (rank == 0 && final_full_matrix != NULL) {
                printf("----------------Final gathered global matrix on root (after sort):----------------\n");
                print_global_matrix_on_root(final_full_matrix, N, rank);
            }
        } else {
            if (rank == 0) {
                fprintf(stderr, "Root: gather_matrix failed for final data (returned %d).\n", gather_res_final);
            }
        }
    }
    
    printf("Sorting completed. Verifying snake-like order...\n");
    int is_sorted = verify_snake_order(&distrib);
    if (rank == 0) {
        if (is_sorted) {
            printf("Verification passed: Matrix is correctly sorted in snake-like order\n");
        } else {
            printf("Verification failed: Matrix is NOT correctly sorted\n");
        }
    }
    
    cleanup_distribution(&distrib);
    if (initial_full_matrix) free(initial_full_matrix);
    if (final_full_matrix) free(final_full_matrix);
    
    MPI_Finalize();
    return 0;
}