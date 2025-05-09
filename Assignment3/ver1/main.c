#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For atoi
#include "quicksort.h" // Includes pivot.h implicitly if quicksort.h does
#include "pivot.h"     // Explicitly include for compare function if not in quicksort.h

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4) {
        if (rank == ROOT) {
            fprintf(stderr, "Usage: %s <input_file> <output_file> <pivot_strategy_id>\n", argv[0]);
            fprintf(stderr, "Pivot strategy IDs:\n");
            fprintf(stderr, "  1: MEDIAN_ROOT\n");
            fprintf(stderr, "  2: MEAN_MEDIAN\n");
            fprintf(stderr, "  3: MEDIAN_MEDIAN\n");
            // Add 0 for SMALLEST_ROOT if you want to expose it, though it's not recommended
        }
        MPI_Finalize();
        return 1;
    }

    char *input_file_name = argv[1];
    char *output_file_name = argv[2];
    int pivot_strategy = atoi(argv[3]);

    if (pivot_strategy < 0 || pivot_strategy > 3) { // Assuming 0 is also a valid internal strategy
         if (rank == ROOT) {
            fprintf(stderr, "Error: Invalid pivot strategy ID %d. Must be 1, 2, or 3.\n", pivot_strategy);
         }
         MPI_Finalize();
         return 1;
    }


    int *all_elements = NULL;
    int n_total = 0; // Total number of elements to sort

    if (rank == ROOT) {
        n_total = read_input(input_file_name, &all_elements);
        if (n_total < 0) {
            fprintf(stderr, "Error reading input file %s on ROOT process.\n", input_file_name);
            // all_elements should be NULL or freed by read_input on error
            MPI_Abort(MPI_COMM_WORLD, 1); // Abort all processes
        }
        if (n_total == 0 && all_elements == NULL) {
             // Handle empty input: write empty output and exit.
            FILE *out_file = fopen(output_file_name, "w");
            if (out_file) {
                // fprintf(out_file, "\n"); // Or just an empty file
                fclose(out_file);
            } else {
                perror("Failed to open output file for empty input");
            }
            printf("0.000000\n"); // Sorting time is 0
            MPI_Finalize();
            return 0;
        }
    }

    // Broadcast total number of elements to all processes
    MPI_Bcast(&n_total, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    if (n_total <= 0 && rank != ROOT) { // Check if other processes should also handle this
        // If ROOT had an issue and aborted, this might not be reached.
        // If n_total is 0 from a valid empty file, other processes proceed with 0 elements.
        if (n_total < 0) { // Should have been caught by MPI_Abort from ROOT
             MPI_Finalize();
             return 1;
        }
    }
    if (n_total == 0) { // All processes handle empty input consistently
        if(rank == ROOT) {
            // Output file already handled by ROOT
            // Time already printed by ROOT
            if(all_elements) free(all_elements); // Should be NULL if n_total is 0
        }
        MPI_Finalize();
        return 0;
    }


    int *my_elements = NULL;
    int my_n; // Number of elements for this process

    // Step 1: Divide the data (distribute from root)
    my_n = distribute_from_root(all_elements, n_total, &my_elements);
    if (my_n < 0) { // Should not happen if distribute_from_root aborts on critical error
        fprintf(stderr, "[Process %d] Error in distributing elements.\n", rank);
        if (my_elements) free(my_elements);
        if (rank == ROOT && all_elements) free(all_elements);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }


    MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes have received data before timing
    double start_time, end_time, sort_time;
    start_time = MPI_Wtime();

    // Step 2: Sort the data locally for each process
    if (my_n > 0 && my_elements != NULL) {
        qsort(my_elements, my_n, sizeof(int), compare);
    }

    // Step 3: Perform global sort
    // The global_sort function modifies my_elements and my_n
    my_n = global_sort(&my_elements, my_n, MPI_COMM_WORLD, pivot_strategy);


    end_time = MPI_Wtime();
    sort_time = end_time - start_time;

    // Gather all elements back on root
    // Note: all_elements on ROOT was used for initial distribution.
    // If n_total is very large, ROOT might need to reallocate all_elements if it was freed,
    // or ensure it's large enough. For now, assume it's still allocated if rank == ROOT.
    // However, it's safer to free it after distribution and re-allocate if needed,
    // or just allocate a new buffer for gathering.
    // For simplicity, let's assume all_elements on ROOT is either still there or needs to be (re)allocated.

    if (rank == ROOT) {
        // If all_elements was freed after distribution, or if its size might have changed
        // due to an uneven initial distribution (though distribute_from_root aims for even),
        // it's safest to ensure it's correctly sized or reallocated.
        // However, gather_on_root expects a buffer of n_total.
        // If all_elements was the original from read_input, it's already n_total.
        // If it was freed, it needs to be reallocated.
        if (!all_elements && n_total > 0) { // If it was freed after distribution
            all_elements = (int *)malloc(n_total * sizeof(int));
            if (!all_elements) {
                fprintf(stderr, "ROOT: Failed to re-allocate all_elements for gather.\n");
                if(my_elements) free(my_elements);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        // If all_elements was potentially modified or is from a different source, ensure it's ready.
        // For this flow, all_elements was from read_input and kept by ROOT.
    }


    gather_on_root(all_elements, my_elements, my_n);

    if (rank == ROOT) {
        // Print sorting time to stdout
        printf("%.6f\n", sort_time);

        // Write sorted data to output file
        if (check_and_print(all_elements, n_total, output_file_name) != 0) {
            fprintf(stderr, "Error writing to output file %s\n", output_file_name);
        }
        if (all_elements) {
            free(all_elements); // Free the initial/gathered array on ROOT
        }
    }

    // Free local elements
    if (my_elements) {
        free(my_elements);
    }

    MPI_Finalize();
    return 0;
}