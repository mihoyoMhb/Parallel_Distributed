// pivot.c
#include "pivot.h"
#include <stdlib.h>
#include <stdio.h>

int compare(const void *v1, const void *v2) {
    return (*(int*)v1 - *(int*)v2);
}

int get_larger_index(int *elements, int n, int val) {
    for (int i = 0; i < n; i++) {
        if (elements[i] > val) return i;
    }
    return n;
}

int get_median(int *elements, int n) {
    return elements[n/2];
}

int select_pivot(int pivot_strategy, int *elements, int n, MPI_Comm communicator) {
    switch(pivot_strategy) {
        case 0: return select_pivot_smallest_root(elements, n, communicator);
        case 1: return select_pivot_median_root(elements, n, communicator);
        case 2: return select_pivot_mean_median(elements, n, communicator);
        case 3: return select_pivot_median_median(elements, n, communicator);
        default: return select_pivot_median_root(elements, n, communicator);
    }
}

int select_pivot_median_root(int *elements, int n, MPI_Comm communicator) {
    int rank, pivot;
    MPI_Comm_rank(communicator, &rank);
    
    if (rank == ROOT) {
        pivot = get_median(elements, n);
    }
    MPI_Bcast(&pivot, 1, MPI_INT, ROOT, communicator);
    return get_larger_index(elements, n, pivot);
}

int select_pivot_mean_median(int *elements, int n, MPI_Comm communicator) {
    int size, rank, med = get_median(elements, n);
    double sum_med;
    
    MPI_Comm_size(communicator, &size);
    MPI_Comm_rank(communicator, &rank);
    
    MPI_Allreduce(&med, &sum_med, 1, MPI_DOUBLE, MPI_SUM, communicator);
    int pivot = (int)(sum_med / size);
    return get_larger_index(elements, n, pivot);
}

int select_pivot_median_median(int *elements, int n, MPI_Comm communicator) {
    int size, rank, med = get_median(elements, n);
    int *all_meds = NULL;
    
    MPI_Comm_size(communicator, &size);
    MPI_Comm_rank(communicator, &rank);
    
    if (rank == ROOT) all_meds = malloc(size * sizeof(int));
    MPI_Gather(&med, 1, MPI_INT, all_meds, 1, MPI_INT, ROOT, communicator);
    
    if (rank == ROOT) {
        qsort(all_meds, size, sizeof(int), compare);
        med = get_median(all_meds, size);
        free(all_meds);
    }
    MPI_Bcast(&med, 1, MPI_INT, ROOT, communicator);
    return get_larger_index(elements, n, med);
}

int select_pivot_smallest_root(int *elements, int n, MPI_Comm communicator) {
    int pivot = elements[0];
    MPI_Bcast(&pivot, 1, MPI_INT, ROOT, communicator);
    return get_larger_index(elements, n, pivot);
}