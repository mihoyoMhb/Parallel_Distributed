#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "data_distribution_v2.h" // Changed from data_distribution.h
#include <errno.h> 

int init_distribution(int N, DataDistribution *distrib) {
    if (distrib == NULL) {
        return -1;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &distrib->world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &distrib->world_rank);
    distrib->N = N;
    distrib->row_comm = MPI_COMM_NULL;
    distrib->col_comm = MPI_COMM_NULL;
    distrib->local_data = NULL;
    return 0;
}

// Modified setup_processor_grid for V2 (1D row decomposition)
int setup_processor_grid(DataDistribution *distrib) {
    int world_size = distrib->world_size;

    // Force p_c = 1, so p_r = world_size (1D decomposition by rows)
    distrib->p_c = 1;
    distrib->p_r = world_size;

    // Calculate this process's position in the grid
    distrib->my_grid_row = distrib->world_rank; // Each process is its own row in the process grid
    distrib->my_grid_col = 0;                   // All processes are in column 0 of the process grid

    // Calculate local block dimensions
    // Each process gets a block of full rows
    if (distrib->p_r == 0) { // Essential check if world_size could be 0, though MPI usually ensures >=1
        if (distrib->world_rank == 0) {
            fprintf(stderr, "Error: Number of processes (p_r) is zero in setup_processor_grid.\n");
        }
        return -1;
    }
    distrib->local_rows = distrib->N / distrib->p_r; 
    distrib->local_cols = distrib->N; // Each process holds full width of its assigned rows

    // Check if N is divisible by p_r (i.e., world_size)
    if (distrib->N % distrib->p_r != 0) {
        if (distrib->world_rank == 0) {
            fprintf(stderr, "Error: Matrix dimension N=%d is not divisible by number of processes p_r=%d.\n",
                   distrib->N, distrib->p_r);
            fprintf(stderr, "This implementation requires N to be a multiple of the total number of processes for 1D row distribution.\n");
        }
        return -1;
    }
    
    if (distrib->local_rows == 0 && distrib->N > 0) { 
        if (distrib->world_rank == 0) {
            fprintf(stderr, "Error: Calculated local_rows is zero (N=%d, p_r=%d). Not enough rows for each process.\n",
                    distrib->N, distrib->p_r);
        }
        return -1; 
    }

    // Create row communicator
    // For p_c=1, my_grid_row = world_rank, my_grid_col = 0.
    // MPI_Comm_split by my_grid_row means each process gets its own row_comm.
    MPI_Comm_split(MPI_COMM_WORLD, distrib->my_grid_row, distrib->my_grid_col, &distrib->row_comm);

    // Create column communicator
    // For p_c=1, my_grid_col = 0 for all.
    // MPI_Comm_split by my_grid_col means all processes form one col_comm.
    MPI_Comm_split(MPI_COMM_WORLD, distrib->my_grid_col, distrib->my_grid_row, &distrib->col_comm);

    if (distrib->world_rank == 0) {
        printf("Created %d x %d process grid (1D row decomposition) for %d x %d matrix.\n", 
               distrib->p_r, distrib->p_c, distrib->N, distrib->N);
        printf("Each process handles %d rows x %d columns (i.e., %d full rows).\n", 
               distrib->local_rows, distrib->local_cols, distrib->local_rows);
    }
    return 0;
}


int allocate_local_block(DataDistribution *distrib) {
    // Ensure local_rows and local_cols are determined and valid from setup_processor_grid
    if (distrib->local_rows < 0 || distrib->local_cols < 0) { // Basic sanity check
        fprintf(stderr, "Process %d: Invalid local dimensions (%d x %d) for allocation.\n",
                distrib->world_rank, distrib->local_rows, distrib->local_cols);
        return -1;
    }

    int local_size = distrib->local_rows * distrib->local_cols;
    
    if (local_size > 0) {
        distrib->local_data = (int*)malloc(local_size * sizeof(int));
        if (distrib->local_data == NULL) {
            fprintf(stderr, "Process %d: Memory allocation failed for local data block of size %d elements.\n",
                    distrib->world_rank, local_size);
            return -1;
        }
    } else {
        distrib->local_data = NULL; // Explicitly NULL for zero size (e.g. N=0 or if local_rows became 0 legitimately)
    }
    return 0;
}

// distribute_matrix should work with 1D decomposition as well.
// Root iterates pr_idx from 0 to p_r-1 (world_size-1)
// pc_idx loop (0 to p_c-1) runs only for pc_idx = 0.
// target_rank = pr_idx * 1 + 0 = pr_idx. So target_rank is essentially the row process index.
// Global indices: global_i = pr_idx * local_rows + i; global_j = 0 * local_cols + j = j;
// This logic needs to be careful for 1D:
// Each process `k` (from 0 to world_size-1) should get rows `k*local_rows` to `(k+1)*local_rows - 1`.
// And it gets all `N` columns for these rows.
int distribute_matrix(int *global_data_on_root, DataDistribution *distrib) {
    int N = distrib->N;
    int world_size = distrib->world_size; // p_r
    int local_rows_per_proc = distrib->local_rows; 
    int cols_per_proc = distrib->local_cols; // Should be N

    if (cols_per_proc != N && N > 0) {
         if(distrib->world_rank == 0) fprintf(stderr, "Distribute_matrix consistency: local_cols (%d) != N (%d)\n", cols_per_proc, N);
         return -1; 
    }

    int *temp_buffer_on_root = NULL;
    int block_size_elements = local_rows_per_proc * cols_per_proc;

    if (distrib->world_rank == 0) {
        if (global_data_on_root != NULL && block_size_elements > 0) {
            temp_buffer_on_root = (int*)malloc(block_size_elements * sizeof(int));
            if (temp_buffer_on_root == NULL) {
                fprintf(stderr, "Root process: Malloc failed for temp_buffer_on_root in distribute_matrix.\n");
                return -1;
            }
        }
    }

    for (int target_proc_rank = 0; target_proc_rank < world_size; target_proc_rank++) {
        if (distrib->world_rank == 0) {
            if (global_data_on_root != NULL) { // Root is sending
                // Pack the block for target_proc_rank
                int global_row_start = target_proc_rank * local_rows_per_proc;
                if (block_size_elements > 0) { // Only pack if there's data
                    int* source_data_ptr;
                    if (target_proc_rank == 0) { // Root preparing its own data
                        // If N=0, block_size_elements is 0, local_data should be NULL or unused.
                        // If N>0, local_data should be allocated.
                        if (!distrib->local_data && block_size_elements > 0) {
                             fprintf(stderr, "Root: local_data is NULL for self-copy, block_size > 0.\n"); return -1;
                        }
                        source_data_ptr = distrib->local_data; 
                        // Pack directly into local_data from global_data_on_root
                        for (int i = 0; i < local_rows_per_proc; i++) {
                            for (int j = 0; j < cols_per_proc; j++) {
                                source_data_ptr[i * cols_per_proc + j] = 
                                    global_data_on_root[(global_row_start + i) * N + j];
                            }
                        }
                    } else { // Root preparing data for another process
                        if (!temp_buffer_on_root && block_size_elements > 0) {
                             fprintf(stderr, "Root: temp_buffer_on_root is NULL for sending to rank %d, block_size > 0.\n", target_proc_rank);
                             // temp_buffer_on_root free is at the end, so just return
                             return -1;
                        }
                        source_data_ptr = temp_buffer_on_root;
                        // Pack into temp_buffer_on_root
                         for (int i = 0; i < local_rows_per_proc; i++) {
                            for (int j = 0; j < cols_per_proc; j++) {
                                source_data_ptr[i * cols_per_proc + j] = 
                                    global_data_on_root[(global_row_start + i) * N + j];
                            }
                        }
                        MPI_Send(source_data_ptr, block_size_elements, MPI_INT, target_proc_rank, 0, MPI_COMM_WORLD);
                    }
                }
            } // else if global_data_on_root is NULL on root, do nothing (e.g. N=0 scenario, or just gathering)
        } else if (distrib->world_rank == target_proc_rank) { // Non-root process is receiving
            if (block_size_elements > 0) { // Only receive if there's data to expect
                if (!distrib->local_data) {
                    fprintf(stderr, "Process %d: local_data is NULL before Recv, expecting %d elements.\n", distrib->world_rank, block_size_elements);
                    return -1; 
                }
                MPI_Recv(distrib->local_data, block_size_elements, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } // else if block_size_elements is 0, do nothing (local_data should be NULL)
        }
    }

    if (distrib->world_rank == 0 && temp_buffer_on_root != NULL) {
        free(temp_buffer_on_root);
    }
    return 0;
}


void local_to_global_indices(int local_row_idx, int local_col_idx,
                           int *global_row_idx, int *global_col_idx,
                           const DataDistribution *distrib) {
    // For 1D row decomposition:
    // my_grid_row is the process rank (0 to p_r-1)
    // local_rows is N / p_r
    // my_grid_col is 0
    // local_cols is N
    *global_row_idx = distrib->my_grid_row * distrib->local_rows + local_row_idx;
    *global_col_idx = local_col_idx; // Since local_cols == N and my_grid_col == 0
}

int global_to_local_indices(int global_row_idx, int global_col_idx,
                          int *local_row_idx, int *local_col_idx,
                          int *owner_proc_grid_row, int *owner_proc_grid_col, // For 1D, owner_proc_grid_col will be 0
                          const DataDistribution *distrib) {
    if (distrib->local_rows == 0 && distrib->N > 0) { // Avoid division by zero if no rows assigned locally, or if N is too small for p_r
        *owner_proc_grid_row = -1; *owner_proc_grid_col = -1;
        *local_row_idx = -1; *local_col_idx = -1;
        return 0; // Cannot belong to this process if local_rows is 0 for this process
    }
    if (distrib->p_r == 0) { // Should be caught by setup, but as a safeguard
        fprintf(stderr, "Error in global_to_local_indices: p_r is zero (rank %d)\n", distrib->world_rank);
        *owner_proc_grid_row = -1; *owner_proc_grid_col = -1;
        *local_row_idx = -1; *local_col_idx = -1;
        return 0;
    }
    
    // For 1D row decomposition (p_c=1, p_r = world_size)
    // local_cols = N
    
    *owner_proc_grid_row = global_row_idx / distrib->local_rows; // This is the rank of the owning process
    *owner_proc_grid_col = 0; // Always 0 for 1D p_c=1 setup

    *local_row_idx = global_row_idx % distrib->local_rows;
    *local_col_idx = global_col_idx; // Global column is local column

    return (*owner_proc_grid_row == distrib->world_rank && *owner_proc_grid_col == distrib->my_grid_col);
}

void cleanup_distribution(DataDistribution *distrib) {
    if (distrib == NULL) {
        return;
    }
    if (distrib->local_data != NULL) {
        free(distrib->local_data);
        distrib->local_data = NULL;
    }
    if (distrib->row_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&distrib->row_comm);
        distrib->row_comm = MPI_COMM_NULL;
    }
    if (distrib->col_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&distrib->col_comm);
        distrib->col_comm = MPI_COMM_NULL;
    }
}

// Corrected read_csv_and_distribute_data function
int read_csv_and_distribute_data(const char* filename, DataDistribution *distrib, int** global_matrix_ptr_for_root_out) {
    int *global_data_on_root = NULL;

    if (distrib->world_rank == 0 && global_matrix_ptr_for_root_out != NULL) {
        *global_matrix_ptr_for_root_out = NULL; 
    }

    if (distrib->world_rank == 0) {
        FILE *fp = fopen(filename, "r");
        if (fp == NULL) {
            fprintf(stderr, "Root: Error opening CSV file %s: %s\n", filename, strerror(errno));
            return -1; 
        }

        int csv_rows, csv_cols;
        if (fscanf(fp, "%d,%d", &csv_rows, &csv_cols) != 2) {
            fprintf(stderr, "Root: Error reading dimensions from CSV file %s\n", filename);
            fclose(fp);
            return -2; 
        }
        // Consume rest of the first line
        int first_line_char; 
        while ((first_line_char = fgetc(fp)) != '\n' && first_line_char != EOF);

        if (csv_rows != distrib->N || csv_cols != distrib->N) {
            fprintf(stderr, "Root: CSV dimensions (%dx%d) do not match N (%dx%d).\n",
                    csv_rows, csv_cols, distrib->N, distrib->N);
            fclose(fp);
            return -3; 
        }

        if (distrib->N > 0) {
            global_data_on_root = (int*)malloc(distrib->N * distrib->N * sizeof(int));
            if (global_data_on_root == NULL) {
                fprintf(stderr, "Root: Failed to allocate memory for global_data_on_root (N=%d)\n", distrib->N);
                fclose(fp);
                return -4; 
            }

            for (int i = 0; i < distrib->N; i++) {
                for (int j = 0; j < distrib->N; j++) {
                    if (fscanf(fp, "%d", &global_data_on_root[i * distrib->N + j]) != 1) {
                        fprintf(stderr, "Root: Error reading matrix data at/before (%d,%d) from CSV\n", i, j);
                        free(global_data_on_root); global_data_on_root = NULL;
                        fclose(fp);
                        return -5; 
                    }
                    // Consume separator or EOL
                    if (j < distrib->N - 1) { // Expect a comma
                        int sep_char = fgetc(fp);
                        while(isspace(sep_char) && sep_char != EOF && sep_char != '\n') sep_char = fgetc(fp); // skip spaces
                        if (sep_char != ',') {
                             fprintf(stderr, "Root: CSV format error after element (%d,%d) - expected ',', got '%c'(%d)\n", i,j, (isprint(sep_char) ? sep_char : '?'), sep_char);
                             free(global_data_on_root); global_data_on_root = NULL;
                             fclose(fp);
                             return -5;
                        }
                    } else { // Expect EOL or EOF after last element of a row
                        int eol_char = fgetc(fp);
                        while(isspace(eol_char) && eol_char != EOF && eol_char != '\n') eol_char = fgetc(fp); // skip spaces
                        if (eol_char != '\n' && eol_char != EOF) {
                            // If it's not newline or EOF (and not just spaces already consumed), it might be an error (e.g. extra data on line)
                            // For simplicity, we can be lenient or strict. Current fscanf will stop at first non-digit for next int.
                            // This simplified version won't strictly check for extra chars after a full row if they are not newlines.
                        }
                         if (eol_char == EOF && i < distrib->N -1) {
                             fprintf(stderr, "Root: Premature EOF in CSV data after row %d\n", i);
                             free(global_data_on_root); global_data_on_root = NULL;
                             fclose(fp);
                             return -5;
                         }
                    }
                }
            }
        } 
        fclose(fp); 

        if (global_matrix_ptr_for_root_out != NULL) {
            *global_matrix_ptr_for_root_out = global_data_on_root;
        }
    }

    int result = distribute_matrix(global_data_on_root, distrib);

    if (distrib->world_rank == 0 && global_data_on_root != NULL) {
        int passed_out = 0;
        if (global_matrix_ptr_for_root_out != NULL && *global_matrix_ptr_for_root_out == global_data_on_root) {
            passed_out = 1;
        }
        if (!passed_out) {
            free(global_data_on_root);
        }
    }

    if (result != 0) {
        return result; 
    }
    
    return 0;
}

// Reverted print_local_block to Ver0 style for concise output
void print_local_block(const DataDistribution *distrib, const char* title) {
    if (!distrib) {
        // Cannot even get rank if distrib is NULL
        printf("Error: print_local_block called with NULL DataDistribution.\n");
        return;
    }

    // Each process prints its own block. Output may interleave without external synchronization.
    printf("Process %d (Grid: %d,%d) - %s (Local Block: %d rows x %d cols, Global N=%d):\n",
           distrib->world_rank, distrib->my_grid_row, distrib->my_grid_col,
           title ? title : "Data",
           distrib->local_rows, distrib->local_cols, distrib->N);

    if (distrib->local_data == NULL) {
        if (distrib->local_rows * distrib->local_cols > 0) {
            printf("  Local data is NULL (but expected %d x %d data).\n", distrib->local_rows, distrib->local_cols);
        } else {
            printf("  Local data is NULL (as expected for 0-size block).\n");
        }
    } else if (distrib->local_rows == 0 || distrib->local_cols == 0) {
        printf("  Local data dimensions are zero or one dimension is zero (local_rows=%d, local_cols=%d).\n",
               distrib->local_rows, distrib->local_cols);
    } else {
        for (int i = 0; i < distrib->local_rows; i++) {
            printf("  Row %d (local): ", i); // Indicate local row index
            for (int j = 0; j < distrib->local_cols; j++) {
                printf("%6d ", distrib->local_data[i * distrib->local_cols + j]);
            }
            printf("\n");
        }
    }
    printf("\n");
    fflush(stdout); // Good practice to flush, especially in parallel programs
}


// gather_full_matrix needs to work correctly with p_c=1 and p_r=world_size
// local_cols = N. local_rows = N/world_size.
// Each process sends its local_rows * N block of data.
// Root needs to place these contiguous blocks one after another.
int gather_full_matrix(DataDistribution *distrib, int *global_matrix_on_root) {
    int N = distrib->N;
    int local_rows = distrib->local_rows; 
    int local_cols = distrib->local_cols; 
    int elements_per_process = local_rows * local_cols;

    if (distrib->world_rank == 0 && N > 0 && global_matrix_on_root == NULL) {
        fprintf(stderr, "Root: gather_full_matrix called with NULL global_matrix_on_root for N=%d.\n", N);
        return -1;
    }
    // elements_per_process < 0 check is removed as local_rows/cols should be >=0 from setup.

    // Ensure local_data is valid or elements_per_process is 0
    if (elements_per_process > 0 && !distrib->local_data) {
        fprintf(stderr, "Process %d: local_data is NULL but trying to send %d elements in gather_full_matrix.\n", distrib->world_rank, elements_per_process);
        return -1; 
    }
    
    // If N=0, elements_per_process will be 0. MPI_Gather with count 0 is valid.
    // local_data can be NULL if elements_per_process is 0.
    
    int mpi_error = MPI_Gather(
        distrib->local_data,        
        elements_per_process,       
        MPI_INT,                    
        global_matrix_on_root,      
        elements_per_process,       
        MPI_INT,                    
        0,                          
        MPI_COMM_WORLD              
    );

    if (mpi_error != MPI_SUCCESS) {
        char err_string[MPI_MAX_ERROR_STRING];
        int len;
        MPI_Error_string(mpi_error, err_string, &len);
        fprintf(stderr, "Process %d: MPI_Gather failed in gather_full_matrix: %s\n", distrib->world_rank, err_string);
        return -2; 
    }

    return 0;
} 