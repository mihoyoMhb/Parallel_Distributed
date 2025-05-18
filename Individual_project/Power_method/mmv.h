#ifndef MMV_H
#define MMV_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// CSR matrix definition
typedef struct {
    int n;          // Global matrix dimension
    int local_n;    // Local number of rows
    int nnz;        // Local non-zero elements
    double *values; // Non-zero element values
    int *col_ind;   // Column indices
    int *row_ptr;   // Row pointers
    int row_start;  // Starting row number for this process
} CSRMatrix;

// Load-balanced distribution based on number of non-zero elements
void distribute_matrix(CSRMatrix *global, CSRMatrix *local, int rank, int size);

// Get the indices of the non-zero elements in the local matrix
/**
 * @brief Identifies and collects unique column indices needed for SpMV computation
 *
 * This function extracts all column indices from a local CSR matrix, sorts them,
 * removes duplicates, and returns the unique set of indices. These indices represent
 * the elements of the x vector that are needed for sparse matrix-vector multiplication.
 *
 * The function performs the following steps:
 * 1. Extracts all column indices from the local matrix
 * 2. Sorts the indices in ascending order using insertion sort
 * 3. Removes duplicate indices
 * 4. Returns the unique set of indices and their count
 *
 * @param local          Pointer to the local CSR matrix containing row pointers, column indices and values
 * @param needed_indices Output parameter that will be allocated and filled with unique column indices
 * @param needed_count   Output parameter that will store the number of unique indices
 *
 * @note The caller is responsible for freeing the memory allocated for needed_indices
 */
void get_needed_indices(CSRMatrix *local, int **needed_indices, int *needed_count);

/**
 * @brief Collects needed x vector elements from all processes and distributes them back
 * 
 * This function implements an optimized communication pattern where:
 * 1. All processes gather their needed indices to process 0
 * 2. Process 0 extracts the corresponding x values from the global vector
 * 3. The needed values are scattered back to each process
 * 
 * @param needed_indices Array of indices needed by this process
 * @param needed_count Number of indices needed by this process
 * @param global_x The global vector x (complete on process 0)
 * @param comm MPI communicator
 * @return Allocated array containing the needed x values for this process
 */
double* collect_and_distribute_x_values(int *needed_indices, int needed_count, 
                                        double *global_x, MPI_Comm comm);

/**
 * @brief Collects local computation results from all processes to form the global result
 * 
 * This function gathers local results from all processes to the root process (rank 0),
 * which assembles them into the complete result vector.
 * 
 * @param local_result The local result vector for this process
 * @param local_size Number of elements in the local result vector
 * @param result Output parameter where the global result is stored (on process 0)
 * @param comm MPI communicator
 */
void collect_results(double *local_result, int local_size, double *result, MPI_Comm comm);

// Optimized SpMV computation - suitable for both single and multi-process scenarios
void spmv_computation(CSRMatrix *local, int *needed_indices, double *needed_x, 
                      int needed_count, double *local_result);

// Optimized serial SpMV computation - using the same optimization strategy as the parallel version
void serial_spmv(CSRMatrix *matrix, double *x, double *result);

// Parallel SpMV computation
void parallel_spmv(CSRMatrix *local, double *global_x, double *result, MPI_Comm comm);

// Load CSR matrix from file
int load_csr_matrix(const char *filename, CSRMatrix *matrix);

#endif 