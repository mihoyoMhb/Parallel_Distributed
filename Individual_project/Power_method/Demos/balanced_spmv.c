#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <time.h>
#include <math.h>

// CSR matrix definition
typedef struct {
    int n;          // Global matrix dimension
    int local_n;    // Local number of rows
    int nnz;        // Local non-zero elements
    double *values; // Non-zero element values
    int *col_ind;   // Column indices
    int *row_ptr;   // Row pointers
    int row_start;  // Starting row number for this process
} CSRMatrix;

int compare_ints(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}
// Load CSR matrix from file
int load_csr_matrix(const char *filename, CSRMatrix *matrix) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        return 0;
    }
    
    int num_rows, num_cols, num_nnz;
    if (fscanf(file, "%d %d %d", &num_rows, &num_cols, &num_nnz) != 3) {
        fclose(file);
        return 0;
    }
    
    matrix->n = num_rows;
    matrix->nnz = num_nnz;
    
    matrix->row_ptr = (int*)malloc((num_rows + 1) * sizeof(int));
    if (num_nnz > 0) {
        matrix->values = (double*)malloc(num_nnz * sizeof(double));
        matrix->col_ind = (int*)malloc(num_nnz * sizeof(int));
    } else {
        matrix->values = NULL;
        matrix->col_ind = NULL;
    }
    
    // Read row_ptr
    for (int i = 0; i <= num_rows; i++) {
        fscanf(file, "%d", &matrix->row_ptr[i]);
    }
    
    // Read col_ind
    for (int i = 0; i < num_nnz; i++) {
        fscanf(file, "%d", &matrix->col_ind[i]);
    }
    
    // Read values
    for (int i = 0; i < num_nnz; i++) {
        fscanf(file, "%lf", &matrix->values[i]);
    }
    
    fclose(file);
    return 1;
}

void distribute_matrix(CSRMatrix *global, CSRMatrix *local, int rank, int size) {
    int *row_starts = NULL;
    
    // Process 0 calculates row distribution
    if (rank == 0) {
        int total_nnz = global->row_ptr[global->n];
        int ideal_distribution = total_nnz / size;
        int remainder = total_nnz % size;
        
        row_starts = (int*)malloc((size + 1) * sizeof(int));
        row_starts[0] = 0;
        
        int current_proc = 0;
        int current_nnz_target = ideal_distribution + (current_proc < remainder ? 1 : 0);
        int nnz_count = 0;
        
        for (int i = 0; i < global->n; i++) {
            // calculate the number of non-zero elements in the current row
            int row_nnz = global->row_ptr[i+1] - global->row_ptr[i];
            // add the number of non-zero elements in the current row to the total number of non-zero elements
            nnz_count += row_nnz;
            
            // if the total number of non-zero elements is greater than or equal to the current target 
            // and the current process is not the last process,
            // increment the current process and set the row start for the current process
            if (nnz_count >= current_nnz_target && current_proc < size - 1) {
                current_proc++;
                row_starts[current_proc] = i + 1;
                current_nnz_target += ideal_distribution + (current_proc < remainder ? 1 : 0);
            }
        }
        row_starts[size] = global->n;
    } else {
        row_starts = (int*)malloc((size + 1) * sizeof(int));
    }
    
    // Broadcast row distribution information to all processes
    MPI_Bcast(row_starts, size + 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Set local matrix parameters
    local->row_start = row_starts[rank];
    // calculate the number of rows for the current process to calculate
    local->local_n = row_starts[rank + 1] - row_starts[rank];
    local->n = global->n;
    local->nnz = global->row_ptr[local->row_start + local->local_n] - global->row_ptr[local->row_start];
    
    // Allocate local matrix memory
    local->values = (double*)malloc(local->nnz * sizeof(double));
    local->col_ind = (int*)malloc(local->nnz * sizeof(int));
    local->row_ptr = (int*)malloc((local->local_n + 1) * sizeof(int));
    
    // Copy data
    if (local->nnz > 0) {
        int start_pos = global->row_ptr[local->row_start];
        memcpy(local->values, &global->values[start_pos], local->nnz * sizeof(double));
        memcpy(local->col_ind, &global->col_ind[start_pos], local->nnz * sizeof(int));
        
        // Set row pointers relative to local matrix
        for (int i = 0; i <= local->local_n; i++) {
            local->row_ptr[i] = global->row_ptr[local->row_start + i] - start_pos;
        }
    } else {
        // If no non-zero elements, all row pointers are 0
        for (int i = 0; i <= local->local_n; i++) {
            local->row_ptr[i] = 0;
        }
    }
    
    free(row_starts);
    
    // Print load balancing information
    printf("Process %d: Processing rows %d to %d (total %d rows, %d non-zero elements)\n", 
           rank, local->row_start, local->row_start + local->local_n - 1,
           local->local_n, local->nnz);
}

// Collect indices of required x vector elements
void get_needed_indices(CSRMatrix *local, int **needed_indices, int *needed_count) {
    // See https://orbit.dtu.dk/files/51272329/tr12_10_Alexandersen_Lazarov_Dammann_1.pdf
    // Page 11-14
    int *temp_indices = (int*)malloc(local->nnz * sizeof(int));
    int count = 0;
    
    // Collect all column indices
    for (int i = 0; i < local->nnz; i++) {
        temp_indices[i] = local->col_ind[i];
    }
    
    // Sort column indices
    // Use qsort to sort the column indices

    qsort(temp_indices, local->nnz, sizeof(int), compare_ints);
    
    // Remove duplicates
    *needed_indices = (int*)malloc(local->nnz * sizeof(int));
    
    if (local->nnz > 0) {
        (*needed_indices)[0] = temp_indices[0];
        count = 1;
        
        for (int i = 1; i < local->nnz; i++) {
            if (temp_indices[i] != temp_indices[i-1]) {
                (*needed_indices)[count] = temp_indices[i];
                count++;
            }
        }
    }
    
    *needed_count = count;
    free(temp_indices);
}

// Optimized SpMV computation - suitable for both single and multi-process scenarios
void spmv_computation(CSRMatrix *local, int *needed_indices, double *needed_x, 
                      int needed_count, double *local_result) {
    // Initialize result vector to 0
    for (int i = 0; i < local->local_n; i++) {
        local_result[i] = 0.0;
    }
    
    // Compute for each row
    for (int i = 0; i < local->local_n; i++) {
        for (int j = local->row_ptr[i]; j < local->row_ptr[i+1]; j++) {
            int col = local->col_ind[j];
            
            // Binary search for the index of the needed x value
            int low = 0, high = needed_count - 1, mid, idx = -1;
            while (low <= high) {
                mid = low + (high - low) / 2;
                if (needed_indices[mid] == col) {
                    idx = mid;
                    break;
                } else if (needed_indices[mid] < col) {
                    low = mid + 1;
                } else {
                    high = mid - 1;
                }
            }
            
            if (idx != -1) {
                local_result[i] += local->values[j] * needed_x[idx];
            }
        }
    }
}

// Optimized serial SpMV computation - using the same optimization strategy as the parallel version
void serial_spmv(CSRMatrix *matrix, double *x, double *result) {
    // Use unified optimization algorithm
    int *needed_indices = NULL;
    int needed_count = 0;
    
    // Find indices of required x elements
    get_needed_indices(matrix, &needed_indices, &needed_count);
    
    // Extract required x element values
    // We don't need every x element, only the ones that are needed for the computation
    double *needed_x = (double*)malloc(needed_count * sizeof(double));
        for (int i = 0; i < needed_count; i++) {
        needed_x[i] = x[needed_indices[i]];
    }
    
    
    // Compute matrix-vector multiplication
    spmv_computation(matrix, needed_indices, needed_x, needed_count, result);
    
    // Free memory
    free(needed_indices);
    free(needed_x);
}

// Parallel SpMV computation
void parallel_spmv(CSRMatrix *local, double *global_x, double *result, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    // If only one process, use optimized serial version
    if (size == 1) {
        serial_spmv(local, global_x, result);
        return;
    }
    
    // Step 1: Determine which x elements each process needs
    int *needed_indices = NULL;
    int needed_count = 0;
    get_needed_indices(local, &needed_indices, &needed_count);
    
    // Step 2: Collect and distribute needed x values
    double *needed_x = (double*)malloc(needed_count * sizeof(double));
    
    // All processes send their needed indices to process 0
    int *recv_counts = NULL;
    int *displs = NULL;
    int *all_indices = NULL;
    
    if (rank == 0) {
        recv_counts = (int*)malloc(size * sizeof(int));
        displs = (int*)malloc(size * sizeof(int));
    }
    
    // Collect how many indices each process needs
    MPI_Gather(&needed_count, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, comm);
    
    // Process 0 computes offsets and allocates space
    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; i++) {
            displs[i] = displs[i-1] + recv_counts[i-1];
        }
        
        int total_indices = displs[size-1] + recv_counts[size-1];
        all_indices = (int*)malloc(total_indices * sizeof(int));
    }
    
    // Collect all indices
    MPI_Gatherv(needed_indices, needed_count, MPI_INT, 
                all_indices, recv_counts, displs, MPI_INT, 0, comm);
    
    // Process 0 prepares corresponding x values
    double *all_x_values = NULL;
    if (rank == 0) {
        int total_indices = displs[size-1] + recv_counts[size-1];
        all_x_values = (double*)malloc(total_indices * sizeof(double));
        
        for (int i = 0; i < total_indices; i++) {
            all_x_values[i] = global_x[all_indices[i]];
        }
    }
    
    // Distribute x values
    MPI_Scatterv(all_x_values, recv_counts, displs, MPI_DOUBLE,
                needed_x, needed_count, MPI_DOUBLE, 0, comm);
    
    // Step 3: Use the unified computation function
    double *local_result = (double*)malloc(local->local_n * sizeof(double));
    spmv_computation(local, needed_indices, needed_x, needed_count, local_result);
    
    // Step 4: Collect results from all processes
    int *result_counts = NULL;
    int *result_displs = NULL;
    
    if (rank == 0) {
        result_counts = (int*)malloc(size * sizeof(int));
        result_displs = (int*)malloc(size * sizeof(int));
    }
    
    // Collect row counts from each process
    MPI_Gather(&local->local_n, 1, MPI_INT, result_counts, 1, MPI_INT, 0, comm);
    
    // Compute result offsets
    if (rank == 0) {
        result_displs[0] = 0;
        for (int i = 1; i < size; i++) {
            result_displs[i] = result_displs[i-1] + result_counts[i-1];
        }
    }
    
    // Collect results
    MPI_Gatherv(local_result, local->local_n, MPI_DOUBLE,
                result, result_counts, result_displs, MPI_DOUBLE, 0, comm);
    
    // Clean up memory
    free(needed_indices);
    free(needed_x);
    free(local_result);
    
    if (rank == 0) {
        free(recv_counts);
        free(displs);
        free(all_indices);
        free(all_x_values);
        free(result_counts);
        free(result_displs);
    }
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    CSRMatrix global_matrix = {0};
    CSRMatrix local_matrix = {0};
    double *x = NULL;
    double *result = NULL;
    double *verification_result = NULL;
    
    // Get matrix filename from command line arguments, or use default
    const char* matrix_filename = (argc > 1) ? argv[1] : "matrix.csr";
    
    // Process 0 reads matrix and generates vector x
    if (rank == 0) {
        if (!loada_csr_matrix(matrix_filename, &global_matrix)) {
            fprintf(stderr, "Error loading matrix from %s\n", matrix_filename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        srand(time(NULL));
        x = (double*)malloc(global_matrix.n * sizeof(double));
        for (int i = 0; i < global_matrix.n; i++) {
            x[i] = (double)rand() / RAND_MAX;
        }
        
        result = (double*)malloc(global_matrix.n * sizeof(double));
        verification_result = (double*)malloc(global_matrix.n * sizeof(double));
    }
    
    // Broadcast matrix dimensions
    MPI_Bcast(&global_matrix.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global_matrix.nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Non-process 0 allocates memory
    if (rank != 0) {
        if (global_matrix.n > 0) {
            global_matrix.row_ptr = (int*)malloc((global_matrix.n + 1) * sizeof(int));
        }
        if (global_matrix.nnz > 0) {
            global_matrix.values = (double*)malloc(global_matrix.nnz * sizeof(double));
            global_matrix.col_ind = (int*)malloc(global_matrix.nnz * sizeof(int));
        }
        
        // Allocate memory for vector x
        x = (double*)malloc(global_matrix.n * sizeof(double));
        result = (double*)malloc(global_matrix.n * sizeof(double));
    }
    
    // Broadcast CSR data
    if (global_matrix.n > 0) {
        MPI_Bcast(global_matrix.row_ptr, global_matrix.n + 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (global_matrix.nnz > 0) {
        MPI_Bcast(global_matrix.values, global_matrix.nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(global_matrix.col_ind, global_matrix.nnz, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    // Broadcast vector x
    MPI_Bcast(x, global_matrix.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Distribute matrix to processes
    distribute_matrix(&global_matrix, &local_matrix, rank, size);
    
    // Compute parallel SpMV
    double start_time = MPI_Wtime();
    parallel_spmv(&local_matrix, x, result, MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    
    // Output and verify results
    if (rank == 0) {
        printf("Parallel computation time: %f seconds\n", end_time - start_time);
        
        printf("Parallel results (first %d elements): ", global_matrix.n < 10 ? global_matrix.n : 10);
        for (int i = 0; i < (global_matrix.n < 10 ? global_matrix.n : 10); i++) {
            printf("%f ", result[i]);
        }
        printf("\n");
        
        // Verify using optimized serial algorithm
        // Prepare global matrix for verification
        CSRMatrix serial_matrix = {0};
        serial_matrix.n = global_matrix.n;
        serial_matrix.local_n = global_matrix.n;
        serial_matrix.nnz = global_matrix.nnz;
        serial_matrix.row_start = 0;
        serial_matrix.values = global_matrix.values;
        serial_matrix.col_ind = global_matrix.col_ind;
        serial_matrix.row_ptr = global_matrix.row_ptr;
        
        // Use the same optimized algorithm for serial computation
        start_time = MPI_Wtime();
        serial_spmv(&serial_matrix, x, verification_result);
        end_time = MPI_Wtime();
        
        printf("Serial optimized computation time: %f seconds\n", end_time - start_time);
        
        printf("Verification results (first %d elements): ", global_matrix.n < 10 ? global_matrix.n : 10);
        for (int i = 0; i < (global_matrix.n < 10 ? global_matrix.n : 10); i++) {
            printf("%f ", verification_result[i]);
        }
        printf("\n");
        
        // Verify consistency between parallel results and optimized serial results
        int errors = 0;
        double tol = 1e-6;
        for (int i = 0; i < global_matrix.n; i++) {
            if (fabs(result[i] - verification_result[i]) > tol) {
                errors++;
                if (errors < 5) {
                    printf("Verification mismatch at index %d: parallel=%f, serial=%f, diff=%e\n",
                            i, result[i], verification_result[i], fabs(result[i] - verification_result[i]));
                }
            }
        }
        
        if (errors > 0) {
            printf("Verification failed with %d mismatches.\n", errors);
        } else {
            printf("Verification passed.\n");
        }
        
        free(verification_result);
    }
    
    // Free memory
    if (local_matrix.values) free(local_matrix.values);
    if (local_matrix.col_ind) free(local_matrix.col_ind);
    if (local_matrix.row_ptr) free(local_matrix.row_ptr);
    
    if (global_matrix.values) free(global_matrix.values);
    if (global_matrix.col_ind) free(global_matrix.col_ind);
    if (global_matrix.row_ptr) free(global_matrix.row_ptr);
    
    if (x) free(x);
    if (result) free(result);
    
    MPI_Finalize();
    return 0;
}