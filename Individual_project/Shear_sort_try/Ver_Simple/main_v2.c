#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <string.h> // For strerror
#include <errno.h>  // For errno

// V2 Includes
#include "data_distribution_v2.h"
#include "row_sort_v2.h"
#include "column_sort_v2.h"
#include "shear_sort_v2.h"

#define DEFAULT_MATRIX_SIZE 16
#define COMPARE_WITH_QSORT 1


// Comparison function for qsort (for standard integer sort)
static int compare_integers(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

// Print a matrix (for debugging purposes)
void print_global_matrix_on_root_v2(int *matrix, int N, int rank) {
    if (rank == 0 && N > 0 && matrix != NULL) {
        printf("Global Matrix (V2 - %dx%d) on Root:\n", N, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%6d ", matrix[i * N + j]);
            }
            printf("\n");
        }
        printf("\n");
    } else if (rank == 0) {
        printf("Global Matrix (V2) on Root: Not printed (N=%d, matrix_ptr=%p)\n", N, (void*)matrix);
    }
}

int main(int argc, char **argv) {
    int N = DEFAULT_MATRIX_SIZE;
    int rank, world_size_main; // Renamed to avoid conflict if DataDistribution has 'size'
    double start_time, end_time;
    const char* csv_filename = "matrix_data.csv"; // Default CSV filename
    int *initial_global_matrix_on_root = NULL; // For rank 0 to get the matrix from read_csv_and_distribute_data

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size_main);

    if (argc > 1) { 
        csv_filename = argv[1]; 
        if (rank == 0) {
            printf("Root: Using CSV file from command line: %s\n", csv_filename);
        }
    } else {
        if (rank == 0) {
            printf("Root: Using default CSV file: %s\n", csv_filename);
        }
    }

    if (rank == 0) {
        FILE *fp = fopen(csv_filename, "r");
        if (fp == NULL) {
            fprintf(stderr, "Root: Error opening CSV file %s: %s. Will attempt to use default N=%d.\n", 
                    csv_filename, strerror(errno), DEFAULT_MATRIX_SIZE);
            // N remains DEFAULT_MATRIX_SIZE
        } else {
            int csv_rows = -1, csv_cols = -1;
            if (fscanf(fp, "%d,%d", &csv_rows, &csv_cols) == 2) {
                if (csv_rows != csv_cols) {
                    fprintf(stderr, "Root: CSV matrix dimensions %dx%d are not square. Using first dimension %d and hoping for the best.\n", 
                            csv_rows, csv_cols, csv_rows);
                }
                if (csv_rows > 0) {
                    N = csv_rows; 
                    printf("Root: Read N=%d from %s\n", N, csv_filename);
                } else {
                    fprintf(stderr, "Root: Invalid dimensions (%d x %d) read from %s. Using default N=%d\n", 
                            csv_rows, csv_cols, csv_filename, DEFAULT_MATRIX_SIZE);
                }
            } else {
                fprintf(stderr, "Root: Error reading dimensions from %s. Using default N=%d\n", 
                        csv_filename, DEFAULT_MATRIX_SIZE);
            }
            fclose(fp);
        }
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (N <= 0) {
        if (rank == 0) fprintf(stderr, "Error: Matrix size N must be positive. N=%d. Aborting.\n", N);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1; 
    }

    DataDistribution distrib;
    if (init_distribution(N, &distrib) != 0) {
        fprintf(stderr, "Process %d: Failed to initialize V2 distribution for N=%d\n", rank, N);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    
    if (setup_processor_grid(&distrib) != 0) { // Calls setup_processor_grid from data_distribution_v2.c
        fprintf(stderr, "Process %d: Failed to setup V2 processor grid\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    
    // Check for invalid local dimensions after setup_processor_grid
    // This check is particularly important for the 1D decomposition
    if (distrib.local_rows <= 0 && N > 0) { // If N > 0, local_rows must be > 0
        if (rank == 0) {
            fprintf(stderr, "Error: Calculated local_rows is %d for N=%d and p_r=%d. Each process must handle at least one row.\n", 
                    distrib.local_rows, N, distrib.p_r);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    if (distrib.local_cols != N && N > 0) { // For 1D row decomp, local_cols must be N
         if (rank == 0) {
            fprintf(stderr, "Error: Consistency check failed. local_cols (%d) != N (%d) for 1D row decomposition.\n", 
                    distrib.local_cols, N);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    if (allocate_local_block(&distrib) != 0) {
        fprintf(stderr, "Process %d: Failed to allocate V2 local block\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    
    int read_dist_result;
    if (rank == 0) {
        read_dist_result = read_csv_and_distribute_data(csv_filename, &distrib, &initial_global_matrix_on_root);
    } else {
        read_dist_result = read_csv_and_distribute_data(csv_filename, &distrib, NULL);
    }
    
    if (read_dist_result != 0) {
        // Error messages are printed inside read_csv_and_distribute_data or distribute_matrix
        if (rank == 0) {
            fprintf(stderr, "Root: Failed to read/distribute data from CSV %s for V2 setup.\n", csv_filename);
            if (initial_global_matrix_on_root != NULL) { // Should be NULL if read_csv_... handled its errors properly
                // free(initial_global_matrix_on_root); // This case indicates an issue in read_csv error path
            }
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }



#if COMPARE_WITH_QSORT
    // Single-process qsort timing on root, using the obtained initial_global_matrix_on_root
    if (rank == 0 && initial_global_matrix_on_root != NULL && N > 0) {
        printf("Root: Performing single-process qsort for timing comparison (on %dx%d matrix)...\n", N, N);
        int* matrix_copy_for_qsort = (int*)malloc((size_t)N * N * sizeof(int));
        if (matrix_copy_for_qsort == NULL) {
            fprintf(stderr, "Root: Failed to allocate memory for qsort copy (N=%d). Skipping qsort timing.\n", N);
        } else {
            memcpy(matrix_copy_for_qsort, initial_global_matrix_on_root, (size_t)N * N * sizeof(int));

            double qsort_start_time = MPI_Wtime();
            qsort(matrix_copy_for_qsort, (size_t)N * N, sizeof(int), compare_integers);
            double qsort_end_time = MPI_Wtime();

            printf("Root: Single-process qsort on %dx%d matrix (standard sort) completed in %f seconds.\n",
                   N, N, qsort_end_time - qsort_start_time);
            
            free(matrix_copy_for_qsort);
        }
        // Free the matrix obtained from read_csv_and_distribute_data, as main is now responsible for it on rank 0
        free(initial_global_matrix_on_root);
        initial_global_matrix_on_root = NULL; 
    } else if (rank == 0 && initial_global_matrix_on_root != NULL && N == 0) {
        // This case should ideally not happen if N=0 implies initial_global_matrix_on_root is NULL from read_csv.
        // But if it somehow is non-NULL, free it.
        free(initial_global_matrix_on_root);
        initial_global_matrix_on_root = NULL;
    }
#endif




    if (N <= 32) { // Print initial local blocks if matrix is small
        print_local_block(&distrib, "Initial local block (V2)");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("Starting Shear Sort V2 (1D Row Decomposition)...\n");
    }

    start_time = MPI_Wtime();
    int pivot_strategy_v2 = 0; // Currently unused in v2 sort functions
    int sort_result = perform_shear_sort_v2(&distrib, pivot_strategy_v2);
    end_time = MPI_Wtime();
    
    if (sort_result != 0) {
        fprintf(stderr, "Process %d: Shear Sort V2 failed\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    
    if (rank == 0) {
        printf("Shear Sort V2 on %dx%d matrix completed in %f seconds\n", N, N, end_time - start_time);
    }

    if (N <= 32) { // Print final local blocks if matrix is small
        print_local_block(&distrib, "Final local block (V2)");
    }
    MPI_Barrier(MPI_COMM_WORLD); 

    int is_sorted_v2 = verify_snake_order_v2(&distrib);
    // verify_snake_order_v2 prints its own PASS/FAIL message on root.
    // It returns 1 if sorted, 0 otherwise (including errors during gather/verify).
    
    cleanup_distribution(&distrib);
    
    MPI_Finalize();
    return (is_sorted_v2 == 1) ? 0 : 1; // Return 0 on success (sorted), 1 on failure (not sorted or error)
} 