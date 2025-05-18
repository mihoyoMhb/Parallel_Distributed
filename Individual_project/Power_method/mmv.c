#include "mmv.h"
#include <math.h>

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
            // https://www.geeksforgeeks.org/binary-search/
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

// Optimized serial SpMV computation
void serial_spmv(CSRMatrix *matrix, double *x, double *result) {
    // Use unified optimization algorithm
    int *needed_indices = NULL;
    int needed_count = 0;
    
    // Find indices of required x elements for non-zero columns
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


double* collect_and_distribute_x_values(int *needed_indices, int needed_count, 
                                        double *global_x, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    // Allocate memory for needed x values
    double *needed_x = (double*)malloc(needed_count * sizeof(double));
    
    // Arrays for communication
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
    
    // Clean up memory
    if (rank == 0) {
        free(recv_counts);
        free(displs);
        free(all_indices);
        free(all_x_values);
    }
    
    return needed_x;
}


void collect_results(double *local_result, int local_size, double *result, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    // Arrays for communication
    int *result_counts = NULL;
    int *result_displs = NULL;
    
    if (rank == 0) {
        result_counts = (int*)malloc(size * sizeof(int));
        result_displs = (int*)malloc(size * sizeof(int));
    }
    
    // Collect row counts from each process
    MPI_Gather(&local_size, 1, MPI_INT, result_counts, 1, MPI_INT, 0, comm);
    
    // Compute result offsets
    if (rank == 0) {
        result_displs[0] = 0;
        for (int i = 1; i < size; i++) {
            result_displs[i] = result_displs[i-1] + result_counts[i-1];
        }
    }
    
    // Collect results
    MPI_Gatherv(local_result, local_size, MPI_DOUBLE,
                result, result_counts, result_displs, MPI_DOUBLE, 0, comm);
    
    // Clean up memory
    if (rank == 0) {
        free(result_counts);
        free(result_displs);
    }
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
    double *needed_x = collect_and_distribute_x_values(needed_indices, needed_count, global_x, comm);
    
    // Step 3: Use the unified computation function
    double *local_result = (double*)malloc(local->local_n * sizeof(double));
    spmv_computation(local, needed_indices, needed_x, needed_count, local_result);
    
    // Step 4: Collect results from all processes
    collect_results(local_result, local->local_n, result, comm);
    
    // Clean up memory
    free(needed_indices);
    free(needed_x);
    free(local_result);
}

// Function to normalize a vector
void normalize_vector(double *vector, int size) {
    double norm = 0.0;
    for (int i = 0; i < size; i++) {
        norm += vector[i] * vector[i];
    }
    norm = sqrt(norm);
    
    if (norm != 0.0) {
        for (int i = 0; i < size; i++) {
            vector[i] /= norm;
        }
    }
}

// Optimized Power method implementation with redundant error checks removed
double power_method(CSRMatrix *matrix, double *initial_vector, int max_iterations, double tolerance, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int n = matrix->n;
    double *current_vector = (double*)malloc(n * sizeof(double));
    double *next_vector = (double*)malloc(n * sizeof(double));

    // Initialize current_vector with the initial_vector or a vector of ones
    if (initial_vector != NULL) {
        memcpy(current_vector, initial_vector, n * sizeof(double));
    } else {
        for (int i = 0; i < n; i++) {
            current_vector[i] = 1.0;
        }
    }


    double eigenvalue = 0.0;
    double prev_eigenvalue = 0.0;
    int converged = 0;

    for (int iter = 0; iter < max_iterations; iter++) {
        // Perform SpMV: next_vector = matrix * current_vector
        if (size > 1) {
            parallel_spmv(matrix, current_vector, next_vector, comm);
        } else {
            serial_spmv(matrix, current_vector, next_vector);
        }

        // Estimate eigenvalue
        eigenvalue = 0.0;
        for (int i = 0; i < n; i++) {
            eigenvalue += current_vector[i] * next_vector[i];
        }
        
        // Aggregate eigenvalue in parallel environment
        if (size > 1) {
            double local_eigenvalue_sum = eigenvalue;
            MPI_Reduce(&local_eigenvalue_sum, &eigenvalue, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            // Broadcast the eigenvalue to all processes
            MPI_Bcast(&eigenvalue, 1, MPI_DOUBLE, 0, comm);
        }

        // Normalize next_vector
        normalize_vector(next_vector, n);

        // Check for convergence (only on rank 0)
        if (rank == 0) {
            converged = (fabs(eigenvalue - prev_eigenvalue) < tolerance) ? 1 : 0;
            
            if (converged) {
                printf("Converged after %d iterations. Eigenvalue: %f\n", iter + 1, eigenvalue);
            } else if (iter == max_iterations - 1) {
                printf("Power method did not converge within %d iterations. Current eigenvalue: %f\n", max_iterations, eigenvalue);
            }
        }
        
        // Broadcast convergence status if parallel
        if (size > 1) {
            MPI_Bcast(&converged, 1, MPI_INT, 0, comm);
        }
            
        if (converged) {
            break;
        }
        
        prev_eigenvalue = eigenvalue;

        // Update current_vector for the next iteration
        memcpy(current_vector, next_vector, n * sizeof(double));
    }

    free(current_vector);
    free(next_vector);

    return eigenvalue;
}



