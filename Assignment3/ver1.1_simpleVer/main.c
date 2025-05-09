#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include "quicksort.h" 
#include "pivot.h"     
// mpiexec -n 1 ./parallel_quicksort /home/mihoyohb/Datas/A3/input500000000.txt /dev/null 2
int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4) {
        if (rank == ROOT) {
            fprintf(stderr, "Usage: %s <input_file> <output_file> <pivot_strategy_id (1, 2, or 3)>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    char *input_file_name = argv[1];
    char *output_file_name = argv[2];
    int pivot_strategy = atoi(argv[3]);

    // Check for pivot strategy parameter validity
    if (pivot_strategy < 1 || pivot_strategy > 3) { 
         if (rank == ROOT) {
            fprintf(stderr, "Error: Invalid pivot strategy ID %d. Must be 1, 2, or 3.\n", pivot_strategy);
         }
         MPI_Finalize();
         return 1;
    }

    int *all_elements = NULL;
    int n_total = 0; 

    if (rank == ROOT) {
        n_total = read_input(input_file_name, &all_elements);
        // read_input handles its own critical memory allocation errors with MPI_Abort
        if (n_total == 0) { 
            FILE *out_file = fopen(output_file_name, "w");
            if (out_file) fclose(out_file); // Assuming success for stable cases
            printf("0.000000\n"); 
            if(all_elements) free(all_elements); 
            MPI_Finalize(); // Finalize ROOT
            // Other processes will get n_total=0 and finalize too.
            // Send a signal or rely on Bcast of n_total=0.
            // Bcast n_total after this block to inform others.
        }
    }

    MPI_Bcast(&n_total, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    if (n_total == 0) { // All processes (including ROOT if it reached here) finalize if n_total is 0.
        MPI_Finalize();
        return 0;
    }

    int *my_elements = NULL;
    int my_n; 

    my_n = distribute_from_root(all_elements, n_total, &my_elements);
    // distribute_from_root handles its own critical memory allocation errors

    MPI_Barrier(MPI_COMM_WORLD); 
    double start_time, end_time, sort_time;
    start_time = MPI_Wtime();

    if (my_n > 0 && my_elements != NULL) {
        qsort(my_elements, my_n, sizeof(int), compare);
    }

    my_n = global_sort(&my_elements, my_n, MPI_COMM_WORLD, pivot_strategy);
    // global_sort handles its own critical memory allocation errors

    end_time = MPI_Wtime();
    sort_time = end_time - start_time;
    
    if (rank == ROOT) {
        // Ensure all_elements is allocated for gathering if it was freed or never allocated (e.g. if n_total was initially 0 then changed)
        // However, if n_total > 0, all_elements was allocated by read_input.
        // If read_input resulted in n_total=0, we would have exited.
        // So, all_elements should be valid if rank == ROOT and n_total > 0.
    }

    gather_on_root(all_elements, my_elements, my_n);
    // gather_on_root handles its own critical memory allocation errors for metadata

    if (rank == ROOT) {
        printf("%.6f\n", sort_time);
        check_and_print(all_elements, n_total, output_file_name);
        // check_and_print assumes file operations succeed
        if (all_elements) {
            free(all_elements); 
        }
    }

    if (my_elements) {
        free(my_elements);
    }

    MPI_Finalize();
    return 0;
}