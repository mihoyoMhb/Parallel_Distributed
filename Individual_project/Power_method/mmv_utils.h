#ifndef MMV_UTILS_H
#define MMV_UTILS_H

#include "mmv.h"
#include <mpi.h>

/**
 * @brief Verifies the results of parallel SpMV computation against serial computation
 * 
 * This function performs the following steps:
 * 1. Prints the parallel computation results
 * 2. Performs a serial computation using the same matrix and vector
 * 3. Prints the serial computation results
 * 4. Compares the parallel and serial results
 * 5. Reports any mismatches
 * 
 * @param global_matrix The global CSR matrix
 * @param x The input vector
 * @param result The result vector from parallel computation
 * @param global_n The global size of the matrix/vectors
 * @return Number of mismatches found (0 if verification passed)
 */
int verify_spmv_results(CSRMatrix *global_matrix, double *x, double *result, int global_n);

#endif /* MMV_UTILS_H */
