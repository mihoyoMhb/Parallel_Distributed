#ifndef ROW_SORT_H
#define ROW_SORT_H

#include <mpi.h>
#include "data_distribution.h"

// Sort direction constants (still used by perform_row_sort logic for choosing compare func)
// Retain if perform_row_sort still conceptually uses them, or if other files might.
// Based on new perform_row_sort, these aren't strictly needed by its signature 
// but the logic determining ascending/descending still exists.
// Let's keep them as they don't harm and reflect sorting intent.
#define SORT_ASCENDING  0
#define SORT_DESCENDING 1

// global_sort declaration is removed as the function is removed from row_sort.c
// int global_sort(int **elements_ptr, int n, MPI_Comm communicator, 
//                int pivot_strategy, int sort_direction);

/**
 * Perform row sorting phase of shearsort
 * 
 * @param distrib         Data distribution structure
 * @param phase_idx       Current phase index (0-based) - (marked as unused in new impl)
 * @param pivot_strategy  Strategy for pivot selection - (marked as unused in new impl)
 * @return                0 on success, non-zero on failure
 */
int perform_row_sort(DataDistribution *distrib, int phase_idx, int pivot_strategy);

#endif /* ROW_SORT_H */