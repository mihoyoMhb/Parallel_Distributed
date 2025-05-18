#ifndef ROW_SORT_V2_H
#define ROW_SORT_V2_H

#include "data_distribution_v2.h"

/**
 * Perform row sort for Shear Sort (V2 - simplified for 1D row decomposition).
 * Each process sorts its assigned full rows locally.
 * Even global rows are sorted ascending, odd global rows descending.
 *
 * @param distrib Pointer to the DataDistribution structure.
 * @param phase_idx The current phase index of the Shear Sort algorithm (determines sort direction if not fixed by global row).
 *                  For Shear Sort, row direction is typically fixed by global row index parity.
 * @return 0 on success, non-zero on error.
 */
int perform_row_sort_v2(DataDistribution *distrib, int phase_idx);

#endif /* ROW_SORT_V2_H */ 