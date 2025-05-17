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
// #define CSV_FILENAME "matrix_data.csv" // Commented out old definition

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
    const char* csv_filename = "matrix_data.csv"; // Default value

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Optional: Allow overriding CSV filename from command line
    if (argc > 1) { 
        csv_filename = argv[1]; 
        if (rank == 0) {
            printf("Root: Using CSV file from command line: %s\n", csv_filename);
        }
    }

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

    int *initial_full_matrix = NULL;
    // This whole block for gathering and printing initial global matrix is conditional on N <= 32
    if (N <= 32) {
        print_local_block(&distrib, "Initial local block");
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize before printing global
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

    int *final_full_matrix = NULL;
    // This whole block for gathering and printing final global matrix is conditional on N <= 32
    
    // printf("Sorting completed. Verifying snake-like order...\n");
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