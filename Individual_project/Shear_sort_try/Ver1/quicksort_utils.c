#include "quicksort_utils.h"

void select_pivot_and_partition(int **elements_ptr, int n, MPI_Comm communicator,
    int *num_small_local, int *num_large_local){
    int rank;
    int k = n;
    MPI_Comm_rank(communicator, &rank);
    int pivot = 0;
    if (rank == 0)
    {
        if(n>0 && *elements_ptr != NULL){
            pivot = (*elements_ptr)[(n - 1)/2];
        }
    }
    for(int i = 0; i < n; i++)
    {
        if((*elements_ptr)[i] > pivot)
        {
            k = i;
            break;
        }
    }
    *num_small_local = k;
    *num_large_local = n - k;
}



int global_sort(int **elements_ptr, int n, MPI_Comm communicator, int is_descending) {
    int my_rank, num_procs;
    MPI_Comm_rank(communicator, &my_rank);
    MPI_Comm_size(communicator, &num_procs);

    int num_small_local, num_large_local;
    select_pivot_and_partition(elements_ptr, n, communicator, &num_small_local, &num_large_local);
    
    // determine process groups and buffers
    
}