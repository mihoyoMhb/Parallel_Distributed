#ifndef COLUMN_SORT_V2_H
#define COLUMN_SORT_V2_H

#include <mpi.h>
#include "data_distribution_v2.h" // Changed include

/**
 * Perform column sorting phase of shearsort (V2)
 * Columns are always sorted in ascending order.
 * This version expects a 1D row decomposition from data_distribution_v2.
 * 
 * @param distrib         Data distribution structure
 * @param phase_idx       Current phase index (0-based) - typically unused for basic column sort
 * @param pivot_strategy  Strategy for pivot selection - typically unused for qsort-based column sort
 * @return                0 on success, non-zero on failure
 */
int perform_column_sort_v2(DataDistribution *distrib, int phase_idx, int pivot_strategy);

#endif /* COLUMN_SORT_V2_H */ 