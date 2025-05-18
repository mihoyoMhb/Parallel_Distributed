#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include "mmv.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    CSRMatrix global_matrix = {0};
    CSRMatrix local_matrix = {0};
    double *x = NULL;
    double *result = NULL;
    double *verification_result = NULL;
    
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
        verification_result = (double*)malloc(global_matrix.n * sizeof(double));
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
    
    // Compute parallel SpMV
    double start_time = MPI_Wtime();
    parallel_spmv(&local_matrix, x, result, MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    
    // Output and verify results
    if (rank == 0) {
        printf("Parallel computation time: %f seconds\n", end_time - start_time);
        
        printf("Parallel results (first %d elements): ", global_matrix.n < 10 ? global_matrix.n : 10);
        for (int i = 0; i < (global_matrix.n < 10 ? global_matrix.n : 10); i++) {
            printf("%f ", result[i]);
        }
        printf("\n");
        
        // Verify using optimized serial algorithm
        // Prepare global matrix for verification
        CSRMatrix serial_matrix = {0};
        serial_matrix.n = global_matrix.n;
        serial_matrix.local_n = global_matrix.n;
        serial_matrix.nnz = global_matrix.nnz;
        serial_matrix.row_start = 0;
        serial_matrix.values = global_matrix.values;
        serial_matrix.col_ind = global_matrix.col_ind;
        serial_matrix.row_ptr = global_matrix.row_ptr;
        
        // Use the same optimized algorithm for serial computation
        start_time = MPI_Wtime();
        serial_spmv(&serial_matrix, x, verification_result);
        end_time = MPI_Wtime();
        
        printf("Serial optimized computation time: %f seconds\n", end_time - start_time);
        
        printf("Verification results (first %d elements): ", global_matrix.n < 10 ? global_matrix.n : 10);
        for (int i = 0; i < (global_matrix.n < 10 ? global_matrix.n : 10); i++) {
            printf("%f ", verification_result[i]);
        }
        printf("\n");
        
        // Verify consistency between parallel results and optimized serial results
        int errors = 0;
        double tol = 1e-6;
        for (int i = 0; i < global_matrix.n; i++) {
            if (fabs(result[i] - verification_result[i]) > tol) {
                errors++;
                if (errors < 5) {
                    printf("Verification mismatch at index %d: parallel=%f, serial=%f, diff=%e\n",
                            i, result[i], verification_result[i], fabs(result[i] - verification_result[i]));
                }
            }
        }
        
        if (errors > 0) {
            printf("Verification failed with %d mismatches.\n", errors);
        } else {
            printf("Verification passed.\n");
        }
        
        free(verification_result);
    }
    
    // Free memory
    if (local_matrix.values) free(local_matrix.values);
    if (local_matrix.col_ind) free(local_matrix.col_ind);
    if (local_matrix.row_ptr) free(local_matrix.row_ptr);
    
    if (global_matrix.values) free(global_matrix.values);
    if (global_matrix.col_ind) free(global_matrix.col_ind);
    if (global_matrix.row_ptr) free(global_matrix.row_ptr);
    
    if (x) free(x);
    if (result) free(result);
    
    MPI_Finalize();
    return 0;
}
