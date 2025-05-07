#include "stencil.h" 
#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define STENCIL_WIDTH 5
#define EXTENT (STENCIL_WIDTH/2)  // radius = 2

int read_input(const char *file_name, double **values);
int write_output(char *file_name, const double *output, int num_values);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s input_file output_file num_steps\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    char *input_file = argv[1];
    char *output_file = argv[2];
    int num_steps = atoi(argv[3]);

    // Only rank 0 reads input data
    double *global_input = NULL;
    int num_values = read_input(input_file, &global_input);


    // Broadcast num_values to all processes
    MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Set up stencil parameters
    double h = 2.0 * PI / num_values;
    const double STENCIL[STENCIL_WIDTH] = {
        1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)
    };

    // Calculate local data size for each process.
    // This one seems useless for requirement in the assignment
    // But we keep it for robustness
    int base = num_values / size;
    int remainder  = num_values % size;
    int local_N = (rank < remainder) ? (base + 1) : base;

    // Set up distribution information for Scatter/Gather (sendcounts and displs)
    int *sendcounts = (int*)malloc(size * sizeof(int));
    int *displacement     = (int*)malloc(size * sizeof(int));
    if (!sendcounts || !displacement) {
        perror("malloc failed");
        return 2;
    }
    int offset = 0;
    for (int i = 0; i < size; i++) {
        sendcounts[i] = (i < remainder) ? (base + 1) : base;
        displacement[i] = offset;
        offset += sendcounts[i];
    }

    // Allocate double buffers, each of size (local_N + 2*EXTENT)
    // Padding is done by adding EXTENT to the left and right of the local data
    // ***Double buffer strategy*** inspired by lecture notes and github code
    // Source: https://github.com/ifromeast/cuda_learning/blob/main/03_gemm/sgemm_v3.cu
    // and https://github.com/xlite-dev/CUDA-Learn-Notes?tab=readme-ov-file
    double *buf_current = (double*)malloc((local_N + 2*EXTENT) * sizeof(double));
    double *buf_next    = (double*)malloc((local_N + 2*EXTENT) * sizeof(double));


    // Distribute global data to each process,
    //    valid data is stored in buf_current[EXTENT, EXTENT+local_N)
    // // Scatter the data to all processes
    // Usage learnt from the tutorial, the link: https://rookiehpc.org/mpi/docs/mpi_scatterv/index.html
    // MPI_Scatterv(a, sendcounts, displacements, MPI_DOUBLE, my_array, my_array_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(
        global_input, sendcounts, displacement, MPI_DOUBLE,
        &buf_current[EXTENT],
        local_N, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );
    if (rank == 0) {
        free(global_input);
        global_input = NULL;
    }

    // Calculate left and right neighbor ranks (periodic boundary)
    int left_rank  = (rank - 1 + size) % size;
    int right_rank = (rank + 1) % size;

    // Define index range of valid data in the buffer
    int eff_start = EXTENT; // valid region starts at EXTENT, we need left 2 elements.
    int eff_end   = EXTENT + local_N;  // valid region length = local_N

    // Ensure timing starts after initial data distribution is complete (timing excludes input and initial distribution)
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // Main loop: perform num_steps iterations of computation
    for (int step = 0; step < num_steps; step++) {
        // Usage: the tutorial link: https://rookiehpc.org/mpi/docs/mpi_irecv/index.html 
        // and https://www.mpich.org/static/docs/v4.1/www3/MPI_Isend.html
        MPI_Request request[4];
        // Receive right pad: store in buf_current[eff_end, eff_end+EXTENT)
        MPI_Irecv(&buf_current[eff_end], EXTENT, MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD, &request[0]);
        // Receive left pad: store in buf_current[0, EXTENT)
        MPI_Irecv(&buf_current[0], EXTENT, MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD, &request[1]);

        // Non-blocking send of own boundary data
        // Send left boundary: buf_current[eff_start, eff_start+EXTENT)
        MPI_Isend(&buf_current[eff_start], EXTENT, MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD, &request[2]);
        // Send right boundary: buf_current[eff_end-EXTENT, eff_end)
        MPI_Isend(&buf_current[eff_end - EXTENT], EXTENT, MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD, &request[3]);

        // Overlap with communication, first compute internal safe region (not dependent on latest padding data),
        // i.e., valid region [eff_start+EXTENT, eff_end-EXTENT)
        for (int i = eff_start + EXTENT; i < eff_end - EXTENT; i++) {
            double sum = 0.0;
            for (int j = -EXTENT; j <= EXTENT; j++) {
                sum += STENCIL[j + EXTENT] * buf_current[i + j];
            }
            buf_next[i] = sum;
        }

        // Wait for non-blocking communication to complete, then compute boundary regions
        MPI_Waitall(4, request, MPI_STATUSES_IGNORE);

        // Compute left boundary region: valid region [eff_start, eff_start+EXTENT)
        // and compute right boundary region: valid region [eff_end-EXTENT, eff_end)
        int offsets[2] = {eff_start, eff_end - EXTENT};
        for (int region = 0; region < 2; region++) {
            for (int idx = 0; idx < EXTENT; idx++) {
                int i = offsets[region] + idx;
                double sum = 0.0;
                for (int j = -EXTENT; j <= EXTENT; j++) {
                    sum += STENCIL[j + EXTENT] * buf_current[i + j];
                }
                buf_next[i] = sum;
            }
        }

        // Using double buffer strategy to swap buf_current and buf_next pointers
        // Origin: lecture notes and provided code
        double *tmp = buf_current;
        buf_current = buf_next;
        buf_next = tmp;
    }

    // After iterations are complete, call Barrier to ensure all processes have finished computation,
    // then stop timing (excluding final data collection time)
    // MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    double local_elapsed = end_time - start_time;
    double max_elapsed;
    // Use MPI_Reduce to find the maximum elapsed time across all processes
    MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Collect final results after timing ends (not included in timing)
    double *final_output = NULL;
    if (rank == 0) {
        final_output = (double*)malloc(num_values * sizeof(double));
        if (!final_output) {
            perror("malloc failed");
            return 2;
        }
    }
    MPI_Gatherv(
        &buf_current[eff_start],  // Only send valid region data
        local_N, MPI_DOUBLE,
        final_output, sendcounts, displacement, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    // Rank 0 outputs timing results, and writes to output file when PRODUCE_OUTPUT_FILE is defined
    if (rank == 0) {
        printf("%f\n", max_elapsed);
#ifdef PRODUCE_OUTPUT_FILE
        if (write_output(output_file, final_output, num_values) != 0) {
            return 2;
        }
#endif
        free(final_output);
    }

    // Free allocated resources
    free(buf_current);
    free(buf_next);
    free(sendcounts);
    free(displacement);

    MPI_Finalize();
    return 0;
}

/* ==========================
 * read_input / write_output
 * ========================== */

int read_input(const char *file_name, double **values) {
    FILE *file;
    if ((file = fopen(file_name, "r")) == NULL) {
        perror("Couldn't open input file");
        return -1;
    }
    int num_values;
    if (fscanf(file, "%d", &num_values) == EOF) {
        perror("Couldn't read element count from input file");
        fclose(file);
        return -1;
    }
    *values = (double*)malloc(num_values * sizeof(double));
    if (!(*values)) {
        perror("malloc for input data failed");
        fclose(file);
        return -1;
    }
    for (int i = 0; i < num_values; i++) {
        if (fscanf(file, "%lf", &((*values)[i])) == EOF) {
            perror("Couldn't read elements from input file");
            fclose(file);
            return -1;
        }
    }
    if (fclose(file) != 0) {
        perror("Warning: couldn't close input file");
    }
    return num_values;
}

int write_output(char *file_name, const double *output, int num_values) {
    FILE *file;
    if ((file = fopen(file_name, "w")) == NULL) {
        perror("Couldn't open output file");
        return -1;
    }
    for (int i = 0; i < num_values; i++) {
        if (fprintf(file, "%.4f ", output[i]) < 0) {
            perror("Couldn't write data");
        }
    }
    if (fprintf(file, "\n") < 0) {
        perror("Couldn't write newline");
    }
    if (fclose(file) != 0) {
        perror("Warning: couldn't close output file");
    }
    return 0;
}