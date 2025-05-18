#include "mmv_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int verify_spmv_results(CSRMatrix *global_matrix, double *x, double *result, int global_n) {
    // Print parallel results
    printf("Parallel results (first %d elements): ", global_n < 10 ? global_n : 10);
    for (int i = 0; i < (global_n < 10 ? global_n : 10); i++) {
        printf("%f ", result[i]);
    }
    printf("\n");
    
    // Verify using optimized serial algorithm
    // Prepare global matrix for verification
    CSRMatrix serial_matrix = {0};
    serial_matrix.n = global_n;
    serial_matrix.local_n = global_n;
    serial_matrix.nnz = global_matrix->nnz;
    serial_matrix.row_start = 0;
    serial_matrix.values = global_matrix->values;
    serial_matrix.col_ind = global_matrix->col_ind;
    serial_matrix.row_ptr = global_matrix->row_ptr;
    
    // Allocate memory for verification result
    double *verification_result = (double*)malloc(global_n * sizeof(double));
    
    // Use the same optimized algorithm for serial computation
    double start_time = MPI_Wtime();
    serial_spmv(&serial_matrix, x, verification_result);
    double end_time = MPI_Wtime();
    
    printf("Serial optimized computation time: %f seconds\n", end_time - start_time);
    
    printf("Verification results (first %d elements): ", global_n < 10 ? global_n : 10);
    for (int i = 0; i < (global_n < 10 ? global_n : 10); i++) {
        printf("%f ", verification_result[i]);
    }
    printf("\n");
    
    // Verify consistency between parallel results and optimized serial results
    int errors = 0;
    double tol = 1e-6;
    for (int i = 0; i < global_n; i++) {
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
    
    return errors;
}
