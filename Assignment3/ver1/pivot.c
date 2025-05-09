#include "pivot.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h> // For qsort, malloc, free, NULL
#include <limits.h> // For INT_MAX if needed

// Compare function for qsort
int compare(const void *v1, const void *v2) {
    int int_v1 = *(const int *)v1;
    int int_v2 = *(const int *)v2;
    if (int_v1 < int_v2) return -1;
    if (int_v1 > int_v2) return 1;
    return 0;
}

// Find the index of the first value in elements that is larger than val
int get_larger_index(int *elements, int n, int val) {
    if (n == 0 || elements == NULL) return 0;
    for (int i = 0; i < n; ++i) {
        if (elements[i] > val) {
            return i;
        }
    }
    return n; // All elements are <= val
}

// Find the median in an array. Assumes array is sorted.
int get_median(int *elements, int n) {
    if (n <= 0 || elements == NULL) {
        // This case should be handled carefully by the calling pivot strategy.
        // Returning 0 might be problematic if 0 is a valid data value.
        // For now, adhering to a simple return, but robust strategies might need more info.
        // fprintf(stderr, "Warning: get_median called with n=%d\n", n);
        return 0; 
    }
    // For even n, this picks the lower of the two middle elements.
    return elements[(n - 1) / 2]; 
}

int select_pivot(int pivot_strategy, int *elements, int n, MPI_Comm communicator) {
    int rank = ROOT; // Default for messages if needed
    MPI_Comm_rank(communicator, &rank);

    switch (pivot_strategy) {
        case 0: // Smallest on root
            return select_pivot_smallest_root(elements, n, communicator);
        case MEDIAN_ROOT: // Defined as 1
            return select_pivot_median_root(elements, n, communicator);
        case MEAN_MEDIAN: // Defined as 2
            return select_pivot_mean_median(elements, n, communicator);
        case MEDIAN_MEDIAN: // Defined as 3
            return select_pivot_median_median(elements, n, communicator);
        default:
            if (rank == ROOT) { // Print error only once
                 fprintf(stderr, "[Process %d] Unknown pivot strategy: %d. Defaulting to MEDIAN_ROOT.\n", rank, pivot_strategy);
            }
            // Fallback to a default. MEDIAN_ROOT is simple.
            return select_pivot_median_root(elements, n, communicator);
    }
}

int select_pivot_median_root(int *elements, int n, MPI_Comm communicator) {
    int rank;
    MPI_Comm_rank(communicator, &rank);

    int pivot_val = 0; // Default pivot value
    if (rank == ROOT) {
        if (n > 0 && elements != NULL) {
            pivot_val = get_median(elements, n); // Assumes elements is sorted
        } else {
            // Root has no elements, pivot_val remains 0 (or some other strategy)
            // fprintf(stderr, "Warning: MEDIAN_ROOT pivot: root (proc 0) has n=%d elements. Pivot is %d.\n", n, pivot_val);
        }
    }

    MPI_Bcast(&pivot_val, 1, MPI_INT, ROOT, communicator);
    return get_larger_index(elements, n, pivot_val);
}

int select_pivot_mean_median(int *elements, int n, MPI_Comm communicator) {
    int rank, size;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    int local_median = 0; 
    if (n > 0 && elements != NULL) {
        local_median = get_median(elements, n); // Assumes elements is sorted
    }

    int *all_medians = NULL;
    if (rank == ROOT) {
        if (size > 0) {
            all_medians = (int *)malloc(size * sizeof(int));
            if (!all_medians) {
                fprintf(stderr, "ROOT: Failed to allocate memory for all_medians in select_pivot_mean_median\n");
                MPI_Abort(communicator, 1);
            }
        }
    }

    MPI_Gather(&local_median, 1, MPI_INT, all_medians, 1, MPI_INT, ROOT, communicator);

    int pivot_val = 0; // Default pivot
    if (rank == ROOT) {
        if (size > 0) {
            long long sum_of_medians = 0; 
            for (int i = 0; i < size; ++i) {
                sum_of_medians += all_medians[i];
            }
            pivot_val = (int)(sum_of_medians / size);
            free(all_medians);
        }
    }

    MPI_Bcast(&pivot_val, 1, MPI_INT, ROOT, communicator);
    return get_larger_index(elements, n, pivot_val);
}

int select_pivot_median_median(int *elements, int n, MPI_Comm communicator) {
    int rank, size;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    int local_median = 0;
    if (n > 0 && elements != NULL) {
        local_median = get_median(elements, n); // Assumes elements is sorted
    }
    
    int *all_medians = NULL;
    if (rank == ROOT) {
        if (size > 0) {
            all_medians = (int *)malloc(size * sizeof(int));
            if (!all_medians) {
                fprintf(stderr, "ROOT: Failed to allocate memory for all_medians in select_pivot_median_median\n");
                MPI_Abort(communicator, 1);
            }
        }
    }

    MPI_Gather(&local_median, 1, MPI_INT, all_medians, 1, MPI_INT, ROOT, communicator);

    int pivot_val = 0; // Default pivot
    if (rank == ROOT) {
        if (size > 0) {
            qsort(all_medians, size, sizeof(int), compare);
            pivot_val = get_median(all_medians, size); // Median of these medians
            free(all_medians);
        }
    }

    MPI_Bcast(&pivot_val, 1, MPI_INT, ROOT, communicator);
    return get_larger_index(elements, n, pivot_val);
}

int select_pivot_smallest_root(int *elements, int n, MPI_Comm communicator) {
    int rank;
    MPI_Comm_rank(communicator, &rank);

    int pivot_val = 0; // Default pivot
    if (rank == ROOT) {
        if (n > 0 && elements != NULL) {
            pivot_val = elements[0]; // Assumes elements is sorted, so smallest is at index 0
        } else {
            // fprintf(stderr, "Warning: SMALLEST_ROOT pivot: root (proc 0) has n=%d. Pivot is %d.\n", n, pivot_val);
        }
    }

    MPI_Bcast(&pivot_val, 1, MPI_INT, ROOT, communicator);
    return get_larger_index(elements, n, pivot_val);
}