#ifndef SHEAR_SORT_H
#define SHEAR_SORT_H

#include <mpi.h>
#include "data_distribution.h"

/**
 * Perform the Shear Sort algorithm on a distributed matrix
 * 
 * @param distrib         Data distribution structure
 * @param pivot_strategy  Strategy for pivot selection
 * @return                0 on success, non-zero on failure
 */
int perform_shear_sort(DataDistribution *distrib, int pivot_strategy);

/**
 * Verify if the matrix is sorted in snake-like order
 * Only root process will have the complete verification result
 * 
 * @param distrib         Data distribution structure
 * @return                1 if correctly sorted, 0 otherwise (on root process)
 */
int verify_snake_order(DataDistribution *distrib);

/**
 * Gather the entire matrix on the root process for verification or output
 * 
 * @param distrib         Data distribution structure
 * @param global_matrix   Buffer to store the gathered matrix (only on root)
 * @return                0 on success, non-zero on failure
 */
int gather_matrix(DataDistribution *distrib, int *global_matrix);

#endif /* SHEAR_SORT_H */