#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include "mmv.h"
#include "mmv_utils.h"

// Configuration parameters
#define VERIFY_SPMV 1  // Set to 0 to disable SpMV verification
#define RUN_POWER_METHOD 1  // Set to 0 to disable power method
#define MAX_ITERATIONS 1000  // Maximum iterations for power method
#define TOLERANCE 1e-6  // Convergence tolerance for power method

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    CSRMatrix global_matrix = {0};
    CSRMatrix local_matrix = {0};
    double *x = NULL;
    double *result = NULL;
    
    // Get matrix filename from command line arguments, or use default
    const char* matrix_filename = (argc > 1) ? argv[1] : "matrix.csr";
    
    // Process 0 reads matrix and generates vector x
    if (rank == 0) {
        if (!load_csr_matrix(matrix_filename, &global_matrix)) {
            fprintf(stderr, "Error loading matrix from %s\n", matrix_filename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        srand(time(NULL));
        x = (double*)malloc(global_matrix.n * sizeof(double));
        for (int i = 0; i < global_matrix.n; i++) {
            x[i] = (double)rand() / RAND_MAX;
        }
        
        result = (double*)malloc(global_matrix.n * sizeof(double));
    }
    
    // Broadcast matrix dimensions
    MPI_Bcast(&global_matrix.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global_matrix.nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Non-process 0 allocates memory
    if (rank != 0) {
        if (global_matrix.n > 0) {
            global_matrix.row_ptr = (int*)malloc((global_matrix.n + 1) * sizeof(int));
        }
        if (global_matrix.nnz > 0) {
            global_matrix.values = (double*)malloc(global_matrix.nnz * sizeof(double));
            global_matrix.col_ind = (int*)malloc(global_matrix.nnz * sizeof(int));
        }
        
        // Allocate memory for vector x
        x = (double*)malloc(global_matrix.n * sizeof(double));
        result = (double*)malloc(global_matrix.n * sizeof(double));
    }
    
    // Broadcast CSR data
    if (global_matrix.n > 0) {
        MPI_Bcast(global_matrix.row_ptr, global_matrix.n + 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (global_matrix.nnz > 0) {
        MPI_Bcast(global_matrix.values, global_matrix.nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(global_matrix.col_ind, global_matrix.nnz, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    // Broadcast vector x
    MPI_Bcast(x, global_matrix.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Distribute matrix to processes
    distribute_matrix(&global_matrix, &local_matrix, rank, size);
    
#if VERIFY_SPMV
    // ===== SpMV verification section =====
    if (rank == 0) {
        printf("\n=== Testing Sparse Matrix-Vector Multiplication ===\n");
    }
    // Normalize the initial vector
    normalize_vector(x, global_matrix.n);
    // Compute parallel SpMV
    double spmv_start_time = MPI_Wtime();
    parallel_spmv(&local_matrix, x, result, MPI_COMM_WORLD);
    double spmv_end_time = MPI_Wtime();
    
    // Output and verify results
    if (rank == 0) {
        printf("Parallel computation time: %f seconds\n", spmv_end_time - spmv_start_time);
        
        // Call the verification function from mmv_utils.c
        verify_spmv_results(&global_matrix, x, result, global_matrix.n);
    }
#endif // VERIFY_SPMV

#if RUN_POWER_METHOD
    // ===== Power Method section =====
    if (rank == 0) {
        printf("\n=== Running Power Method for Dominant Eigenvector ===\n");
        
        // Reset the x vector for power method
        for (int i = 0; i < global_matrix.n; i++) {
            x[i] = (double)rand() / RAND_MAX;
        }
    }
    
    // Broadcast the initial vector
    MPI_Bcast(x, global_matrix.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Run power method and measure time
    double power_start_time = MPI_Wtime();
    double eigenvalue = power_method(&local_matrix, x, MAX_ITERATIONS, TOLERANCE, MPI_COMM_WORLD);
    double power_end_time = MPI_Wtime();
    
    if (rank == 0) {
        printf("Eigenvalue: %f\n", eigenvalue);
        printf("Power method computation time: %f seconds\n", power_end_time - power_start_time);
        // Print a sample of the eigenvector
        printf("Dominant eigenvector (first %d elements): ", global_matrix.n < 10 ? global_matrix.n : 10);
        for (int i = 0; i < (global_matrix.n < 10 ? global_matrix.n : 10); i++) {
            printf("%f ", x[i]);
        }
        printf("\n");
    }
#endif // RUN_POWER_METHOD
    
    // Free resources
    if (x) free(x);
    if (result) free(result);
    
    // Free matrix memory
    if (global_matrix.row_ptr) free(global_matrix.row_ptr);
    if (global_matrix.col_ind) free(global_matrix.col_ind);
    if (global_matrix.values) free(global_matrix.values);
    
    if (local_matrix.row_ptr) free(local_matrix.row_ptr);
    if (local_matrix.col_ind) free(local_matrix.col_ind);
    if (local_matrix.values) free(local_matrix.values);
    
    MPI_Finalize();
    return 0;
}
