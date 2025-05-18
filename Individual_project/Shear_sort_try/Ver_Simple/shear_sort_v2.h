#ifndef SHEAR_SORT_V2_H
#define SHEAR_SORT_V2_H

#include <mpi.h>
#include "data_distribution_v2.h"

/**
 * Perform the Shear Sort algorithm (V2) on a distributed matrix.
 * Uses simplified row sort and standard column sort, expecting 1D row decomposition.
 * 
 * @param distrib         Data distribution structure (from data_distribution_v2.h)
 * @param pivot_strategy  (Currently unused in v2 sort functions, but kept for signature compatibility)
 * @return                0 on success, non-zero on failure
 */
int perform_shear_sort_v2(DataDistribution *distrib, int pivot_strategy);

/**
 * Verify if the matrix is sorted in snake-like order (V2).
 * This function will gather the matrix to the root process for verification.
 * Uses the gather_matrix_v2 function internally.
 * Only root process will have the complete verification result.
 * 
 * @param distrib         Data distribution structure (from data_distribution_v2.h)
 * @return                1 if correctly sorted, 0 otherwise (on root process); -1 on error during gather.
 */
int verify_snake_order_v2(DataDistribution *distrib);

// Note: The gather_matrix function is now part of data_distribution_v2.c as gather_full_matrix.
// verify_snake_order_v2 will use that.

#endif /* SHEAR_SORT_V2_H */ 