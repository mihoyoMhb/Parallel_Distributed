#ifndef COLUMN_SORT_H
#define COLUMN_SORT_H

#include <mpi.h>
#include "data_distribution.h"

/**
 * Perform column sorting phase of shearsort
 * 
 * @param distrib         Data distribution structure
 * @param phase_idx       Current phase index (0-based)
 * @param pivot_strategy  Strategy for pivot selection
 * @return                0 on success, non-zero on failure
 */
int perform_column_sort(DataDistribution *distrib, int phase_idx, int pivot_strategy);

#endif /* COLUMN_SORT_H */