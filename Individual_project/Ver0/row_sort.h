#ifndef ROW_SORT_H
#define ROW_SORT_H

#include <mpi.h>
#include "data_distribution.h"

// Sort direction constants
#define SORT_ASCENDING  0
#define SORT_DESCENDING 1

/**
 * Perform global parallel sorting on a distributed array segment
 * 
 * @param elements_ptr    Pointer to elements array (may be reallocated)
 * @param n               Number of elements
 * @param communicator    MPI communicator for the sort operation
 * @param pivot_strategy  Strategy for pivot selection
 * @param sort_direction  Sort direction (SORT_ASCENDING or SORT_DESCENDING)
 * @return               Number of elements after sorting
 */
int global_sort(int **elements_ptr, int n, MPI_Comm communicator, 
               int pivot_strategy, int sort_direction);

/**
 * Perform row sorting phase of shearsort
 * 
 * @param distrib         Data distribution structure
 * @param phase_idx       Current phase index (0-based)
 * @param pivot_strategy  Strategy for pivot selection
 * @return                0 on success, non-zero on failure
 */
int perform_row_sort(DataDistribution *distrib, int phase_idx, int pivot_strategy);

#endif /* ROW_SORT_H */