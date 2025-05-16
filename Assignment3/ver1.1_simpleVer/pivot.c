#include "pivot.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h> 
#include <limits.h> 

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
    return n; 
}

// Find the median in an array. Assumes array is sorted.
int get_median(int *elements, int n) {
    if (n <= 0 || elements == NULL) {
        return 0; 
    }
    return elements[(n - 1) / 2]; 
}

int select_pivot(int pivot_strategy, int *elements, int n, MPI_Comm communicator) {
    switch (pivot_strategy) {
        case 0: 
            return select_pivot_smallest_root(elements, n, communicator);
        case MEDIAN_ROOT: 
            return select_pivot_median_root(elements, n, communicator);
        case MEAN_MEDIAN: 
            return select_pivot_mean_median(elements, n, communicator);
        case MEDIAN_MEDIAN: 
            return select_pivot_median_median(elements, n, communicator);
        default:
            return select_pivot_median_root(elements, n, communicator);
    }
}

int select_pivot_median_root(int *elements, int n, MPI_Comm communicator) {
    int rank;
    MPI_Comm_rank(communicator, &rank);
    int pivot_val = 0; 
    if (rank == ROOT) {
        if (n > 0 && elements != NULL) {
            pivot_val = get_median(elements, n); 
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
        local_median = get_median(elements, n); 
    }

    int *all_medians = NULL;
    if (rank == ROOT && size > 0) { 
        all_medians = (int *)malloc(size * sizeof(int));
    }

    MPI_Gather(&local_median, 1, MPI_INT, all_medians, 1, MPI_INT, ROOT, communicator);

    int pivot_val = 0; 
    if (rank == ROOT) {
        long long sum_of_medians = 0;
        for (int i = 0; i < size; ++i) {
            sum_of_medians += all_medians[i];
        }
        pivot_val = (int)(sum_of_medians / size);
        free(all_medians);
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
        local_median = get_median(elements, n); 
    }
    
    int *all_medians = NULL;
    if (rank == ROOT && size > 0) { // Only allocate if communicator is not empty
        all_medians = (int *)malloc(size * sizeof(int));
    }

    MPI_Gather(&local_median, 1, MPI_INT, all_medians, 1, MPI_INT, ROOT, communicator);

    int pivot_val = 0; 
    if (rank == ROOT) {
        qsort(all_medians, size, sizeof(int), compare);
        pivot_val = get_median(all_medians, size); 
        free(all_medians);
    }

    MPI_Bcast(&pivot_val, 1, MPI_INT, ROOT, communicator);
    return get_larger_index(elements, n, pivot_val);
}

int select_pivot_smallest_root(int *elements, int n, MPI_Comm communicator) {
    return 0; // This function is not implemented in the original code
}