#include "stencil.h"
#include <mpi.h>
#include <string.h>

int main(int argc, char **argv) {
    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (4 != argc) {
        if (world_rank == 0) {
            printf("Usage: stencil input_file output_file number_of_applications\n");
        }
        MPI_Finalize();
        return 1;
    }
    char *input_name = argv[1];
    char *output_name = argv[2];
    int num_steps = atoi(argv[3]);

    // Only rank 0 reads the input file
    double *global_input = NULL;
    int num_values = 0;
    if (world_rank == 0) {
        if (0 > (num_values = read_input(input_name, &global_input))) {
            MPI_Abort(MPI_COMM_WORLD, 2);
            return 2;
        }
    }

    // Broadcast number of values to all processes
    MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Stencil values
    double h = 2.0*PI/num_values;
    const int STENCIL_WIDTH = 5;
    const int EXTENT = STENCIL_WIDTH/2;
    const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};

    // Calculate local data size for each process
    int base = num_values / world_size;
    int remainder = num_values % world_size;
    int local_N = (world_rank < remainder) ? base + 1 : base;

    // Calculate displacements and counts for scatter/gather operations
    int *sendcounts = malloc(world_size * sizeof(int));
    int *displs = malloc(world_size * sizeof(int));
    
    int offset = 0;
    for (int i = 0; i < world_size; i++) {
        sendcounts[i] = (i < remainder) ? base + 1 : base;
        displs[i] = offset;
        offset += sendcounts[i];
    }

    // Allocate memory for local input and output arrays
    double *local_input = malloc(local_N * sizeof(double));
    double *local_output = malloc(local_N * sizeof(double));
    
    // Temporary arrays for full data (each process will have a copy)
    double *input = malloc(num_values * sizeof(double));
    double *output = malloc(num_values * sizeof(double));
    
    if (!local_input || !local_output || !input || !output) {
        perror("Memory allocation failed");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    // FIXED: First copy data to input on rank 0, then broadcast
    if (world_rank == 0) {
        memcpy(input, global_input, num_values * sizeof(double));
        free(global_input); // Free after copy, no longer needed
    }
    
    // Now broadcast the initialized input data to all processes
    MPI_Bcast(input, num_values, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Start timer
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    // Repeatedly apply stencil
    for (int s = 0; s < num_steps; s++) {
        // Distribute work among processes:
        // Each process computes a portion of the array
        int start_idx = displs[world_rank];
        int end_idx = start_idx + local_N;
        
        // Apply stencil to local section
        for (int i = start_idx; i < end_idx; i++) {
            double result = 0;
            
            // Use the appropriate calculation method based on the index
            if (i < EXTENT) {
                // Left boundary with periodic wrapping
                for (int j = 0; j < STENCIL_WIDTH; j++) {
                    int index = (i - EXTENT + j + num_values) % num_values;
                    result += STENCIL[j] * input[index];
                }
            } 
            else if (i >= num_values - EXTENT) {
                // Right boundary with periodic wrapping
                for (int j = 0; j < STENCIL_WIDTH; j++) {
                    int index = (i - EXTENT + j) % num_values;
                    result += STENCIL[j] * input[index];
                }
            } 
            else {
                // Interior points
                for (int j = 0; j < STENCIL_WIDTH; j++) {
                    int index = i - EXTENT + j;
                    result += STENCIL[j] * input[index];
                }
            }
            
            // Store the result in local output
            local_output[i - start_idx] = result;
        }
        
        // Gather all local results to create the complete output array
        MPI_Allgatherv(local_output, local_N, MPI_DOUBLE, 
                      output, sendcounts, displs, MPI_DOUBLE, 
                      MPI_COMM_WORLD);
        
        // Swap input and output for next iteration
        if (s < num_steps - 1) {
            double *tmp = input;
            input = output;
            output = tmp;
        }
    }

    // Stop timer
    double my_execution_time = MPI_Wtime() - start;
    double max_time;
    MPI_Reduce(&my_execution_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    // Only rank 0 writes output and prints time
    if (world_rank == 0) {
        printf("%f\n", max_time);
#ifdef PRODUCE_OUTPUT_FILE
        if (0 != write_output(output_name, output, num_values)) {
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
#endif
    }
    
    // Clean up
    free(local_input);
    free(local_output);
    free(input);
    free(output);
    free(sendcounts);
    free(displs);
    
    MPI_Finalize();
    return 0;
}

int read_input(const char *file_name, double **values) {
    FILE *file;
    if (NULL == (file = fopen(file_name, "r"))) {
        perror("Couldn't open input file");
        return -1;
    }
    int num_values;
    if (EOF == fscanf(file, "%d", &num_values)) {
        perror("Couldn't read element count from input file");
        return -1;
    }
    if (NULL == (*values = malloc(num_values * sizeof(double)))) {
        perror("Couldn't allocate memory for input");
        return -1;
    }
    for (int i=0; i<num_values; i++) {
        if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
            perror("Couldn't read elements from input file");
            return -1;
        }
    }
    if (0 != fclose(file)){
        perror("Warning: couldn't close input file");
    }
    return num_values;
}

int write_output(char *file_name, const double *output, int num_values) {
    FILE *file;
    if (NULL == (file = fopen(file_name, "w"))) {
        perror("Couldn't open output file");
        return -1;
    }
    for (int i = 0; i < num_values; i++) {
        if (0 > fprintf(file, "%.4f ", output[i])) {
            perror("Couldn't write to output file");
        }
    }
    if (0 > fprintf(file, "\n")) {
        perror("Couldn't write to output file");
    }
    if (0 != fclose(file)) {
        perror("Warning: couldn't close output file");
    }
    return 0;
}