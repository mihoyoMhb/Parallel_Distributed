#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "data_distribution.h"
#include <errno.h> // For errno

int init_distribution(int N, DataDistribution *distrib) {
    if (distrib == NULL) {
        return -1;
    }
    
    // Initialize basic MPI information
    MPI_Comm_size(MPI_COMM_WORLD, &distrib->world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &distrib->world_rank);
    
    // Store global matrix size
    distrib->N = N;
    
    // Initialize communicator variables
    distrib->row_comm = MPI_COMM_NULL;
    distrib->col_comm = MPI_COMM_NULL;
    
    // Initialize data pointer to NULL
    distrib->local_data = NULL;
    
    return 0;
}

int setup_processor_grid(DataDistribution *distrib) {
    // Determine grid dimensions based on total processes
    // Try to make grid as square as possible
    int world_size = distrib->world_size;
    
    // Find factors of world_size that are as close as possible to sqrt(world_size)
    int p_r = (int)sqrt(world_size);
    while (world_size % p_r != 0 && p_r > 1) {
        p_r--;
    }
    int p_c = world_size / p_r;
    
    distrib->p_r = p_r;
    distrib->p_c = p_c;
    
    // Calculate this process's position in the grid
    distrib->my_grid_row = distrib->world_rank / p_c;
    distrib->my_grid_col = distrib->world_rank % p_c;
    
    // Calculate local block dimensions
    distrib->local_rows = distrib->N / p_r;
    distrib->local_cols = distrib->N / p_c;
    
    // Check if N is divisible by p_r and p_c
    if (distrib->N % p_r != 0 || distrib->N % p_c != 0) {
        if (distrib->world_rank == 0) {
            printf("Warning: Matrix dimension N=%d is not divisible by grid dimensions (%d,%d)\n", 
                   distrib->N, p_r, p_c);
            printf("This implementation requires even distribution of data.\n");
        }
        return -1;
    }
    
    // Create row communicator (processes with the same row coordinate)
    MPI_Comm_split(MPI_COMM_WORLD, distrib->my_grid_row, distrib->my_grid_col, &distrib->row_comm);
    
    // Create column communicator (processes with the same column coordinate)
    MPI_Comm_split(MPI_COMM_WORLD, distrib->my_grid_col, distrib->my_grid_row, &distrib->col_comm);
    
    if (distrib->world_rank == 0) {
        printf("Created %d×%d process grid for %d×%d matrix.\n", p_r, p_c, distrib->N, distrib->N);
        printf("Each process handles %d×%d local block.\n", distrib->local_rows, distrib->local_cols);
    }
    
    return 0;
}

int allocate_local_block(DataDistribution *distrib) {
    int local_size = distrib->local_rows * distrib->local_cols;
    
    // Allocate memory for local data block
    distrib->local_data = (int*)malloc(local_size * sizeof(int));
    if (distrib->local_data == NULL) {
        fprintf(stderr, "Process %d: Memory allocation failed for local data block\n", 
                distrib->world_rank);
        return -1;
    }
    
    return 0;
}

int distribute_matrix(int *global_data, DataDistribution *distrib) {
    int N = distrib->N;
    int p_r = distrib->p_r;
    int p_c = distrib->p_c;
    int local_rows = distrib->local_rows;
    int local_cols = distrib->local_cols;
    
    // Temporary buffer for sending/receiving blocks
    int *temp_buffer = NULL;
    
    if (distrib->world_rank == 0) {
        temp_buffer = (int*)malloc(local_rows * local_cols * sizeof(int));
        if (temp_buffer == NULL) {
            fprintf(stderr, "Root process: Memory allocation failed for temp buffer\n");
            return -1;
        }
    }
    
    // Root process distributes data to all processes
    for (int pr = 0; pr < p_r; pr++) {
        for (int pc = 0; pc < p_c; pc++) {
            int target_rank = pr * p_c + pc;
            
            // If this is the target process, just copy data locally
            if (target_rank == 0 && distrib->world_rank == 0) {
                for (int i = 0; i < local_rows; i++) {
                    for (int j = 0; j < local_cols; j++) {
                        int global_i = pr * local_rows + i;
                        int global_j = pc * local_cols + j;
                        distrib->local_data[i * local_cols + j] = global_data[global_i * N + global_j];
                    }
                }
            }
            // If this is the root process, send to other processes
            else if (distrib->world_rank == 0) {
                // Pack data for target process
                for (int i = 0; i < local_rows; i++) {
                    for (int j = 0; j < local_cols; j++) {
                        int global_i = pr * local_rows + i;
                        int global_j = pc * local_cols + j;
                        temp_buffer[i * local_cols + j] = global_data[global_i * N + global_j];
                    }
                }
                
                // Send the packed data to target process
                MPI_Send(temp_buffer, local_rows * local_cols, MPI_INT, target_rank, 0, MPI_COMM_WORLD);
            }
            // If this is a receiving process
            else if (distrib->world_rank == target_rank) {
                // Receive data from root
                MPI_Recv(distrib->local_data, local_rows * local_cols, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    
    // Free temporary buffer
    if (distrib->world_rank == 0 && temp_buffer != NULL) {
        free(temp_buffer);
    }
    
    return 0;
}

void local_to_global_indices(int local_row, int local_col, 
                           int *global_row, int *global_col, 
                           const DataDistribution *distrib) {
    *global_row = distrib->my_grid_row * distrib->local_rows + local_row;
    *global_col = distrib->my_grid_col * distrib->local_cols + local_col;
}

int global_to_local_indices(int global_row, int global_col,
                          int *local_row, int *local_col,
                          int *proc_row, int *proc_col,
                          const DataDistribution *distrib) {
    // Calculate which process owns this global index
    *proc_row = global_row / distrib->local_rows;
    *proc_col = global_col / distrib->local_cols;
    
    // Calculate the local indices within that process
    *local_row = global_row % distrib->local_rows;
    *local_col = global_col % distrib->local_cols;
    
    // Check if this index belongs to the current process
    return (*proc_row == distrib->my_grid_row && *proc_col == distrib->my_grid_col);
}

void cleanup_distribution(DataDistribution *distrib) {
    if (distrib == NULL) {
        return;
    }
    
    // Free local data
    if (distrib->local_data != NULL) {
        free(distrib->local_data);
        distrib->local_data = NULL;
    }
    
    // Free communicators
    if (distrib->row_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&distrib->row_comm);
    }
    
    if (distrib->col_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&distrib->col_comm);
    }
}

int read_csv_and_distribute_data(const char* filename, DataDistribution *distrib) {
    if (distrib->world_rank == 0) {
        FILE *fp = fopen(filename, "r");
        if (fp == NULL) {
            fprintf(stderr, "Root: Error opening CSV file %s: %s\n", filename, strerror(errno));
            return -1;
        }

        int csv_rows, csv_cols;
        // Read dimensions from the first line
        if (fscanf(fp, "%d,%d\n", &csv_rows, &csv_cols) != 2) {
            fprintf(stderr, "Root: Error reading dimensions from CSV file %s\n", filename);
            fclose(fp);
            return -2;
        }

        // Basic check: Ensure CSV dimensions match the initialized N
        // More robustly, N should be determined FROM the CSV by root and broadcasted before init_distribution
        // For this implementation, we assume N was correctly set based on CSV prior to this call.
        if (csv_rows != distrib->N || csv_cols != distrib->N) {
            fprintf(stderr, "Root: CSV dimensions (%d x %d) do not match expected N (%d x %d)\n", 
                    csv_rows, csv_cols, distrib->N, distrib->N);
            fclose(fp);
            return -3; // Or handle resizing/reinitialization if N was a default
        }

        int *global_data = (int*)malloc(distrib->N * distrib->N * sizeof(int));
        if (global_data == NULL) {
            fprintf(stderr, "Root: Failed to allocate memory for global_data to read CSV\n");
            fclose(fp);
            return -4;
        }

        for (int i = 0; i < distrib->N; i++) {
            for (int j = 0; j < distrib->N; j++) {
                if (fscanf(fp, "%d", &global_data[i * distrib->N + j]) != 1) {
                    fprintf(stderr, "Root: Error reading matrix data at (%d,%d) from CSV\n", i, j);
                    free(global_data);
                    fclose(fp);
                    return -5;
                }
                // Consume comma or newline
                if (j < distrib->N - 1) {
                    if (fgetc(fp) != ',') {
                        // allow for optional trailing comma on line or space
                    }
                } else {
                    while (fgetc(fp) != '\n' && !feof(fp)); // consume rest of line
                }
            }
        }
        fclose(fp);

        // Now distribute this global_data using the existing distribute_matrix function
        int result = distribute_matrix(global_data, distrib);
        free(global_data);
        if (result != 0) {
             fprintf(stderr, "Root: Failed during distribute_matrix after CSV read\n");
            return -6;
        }
        printf("Root process successfully read %s and initiated data distribution.\n", filename);

    } else {
        // Non-root processes will receive their data via distribute_matrix, called by root.
        // So they just need to participate in the collective operations within distribute_matrix.
        // This means calling distribute_matrix with NULL for global_data.
        int result = distribute_matrix(NULL, distrib); // Matched with root's call
        if (result != 0) {
            fprintf(stderr, "Process %d: Failed during distribute_matrix (receiving part)\n", distrib->world_rank);
            return -7;
        }
    }
    return 0;
}

void print_local_block(const DataDistribution *distrib, const char* title) {
    printf("Process %d (%d,%d) - %s (%d x %d local block):\n", 
           distrib->world_rank, distrib->my_grid_row, distrib->my_grid_col, 
           title, distrib->local_rows, distrib->local_cols);
    if (distrib->local_data == NULL) {
        printf("  Local data is NULL.\n");
        return;
    }
    for (int i = 0; i < distrib->local_rows; i++) {
        printf("  ");
        for (int j = 0; j < distrib->local_cols; j++) {
            printf("%6d ", distrib->local_data[i * distrib->local_cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int gather_full_matrix(DataDistribution *distrib, int *global_matrix) {
    int N = distrib->N;
    int p_r = distrib->p_r;
    int p_c = distrib->p_c;
    int local_rows = distrib->local_rows;
    int local_cols = distrib->local_cols;
    int local_size = local_rows * local_cols;

    if (distrib->world_rank == 0 && global_matrix == NULL) {
        fprintf(stderr, "Root: gather_full_matrix called with NULL global_matrix buffer.\n");
        return -1;
    }

    // Each process sends its local_data to the root.
    // Root receives these blocks and assembles them.

    if (distrib->world_rank == 0) {
        // Root process itself copies its own data first
        for (int i = 0; i < local_rows; i++) {
            for (int j = 0; j < local_cols; j++) {
                int global_i = distrib->my_grid_row * local_rows + i; // Should be 0 for root
                int global_j = distrib->my_grid_col * local_cols + j; // Should be 0 for root
                global_matrix[global_i * N + global_j] = distrib->local_data[i * local_cols + j];
            }
        }

        // Receive data from other processes
        int *temp_buffer = (int*)malloc(local_size * sizeof(int));
        if (temp_buffer == NULL && local_size > 0) { // only error if malloc fails for non-empty buffer
            fprintf(stderr, "Root: Malloc failed for temp_buffer in gather_full_matrix\n");
            // Continue without it for rank 0? Or abort? For now, it might only affect receiving.
            // If world_size is 1, this buffer is not used for receiving.
        }

        for (int r = 0; r < p_r; r++) {
            for (int c = 0; c < p_c; c++) {
                int source_rank = r * p_c + c;
                if (source_rank == 0) continue; // Already handled root's own data

                // MPI_Status status;
                int recv_count = 0;
                // Use MPI_Probe to find out the actual size of the message from this specific source
                // This is not strictly necessary if all blocks are guaranteed to be local_rows * local_cols
                // But it's safer if block sizes could vary (they don't in current setup)
                // MPI_Probe(source_rank, 0, MPI_COMM_WORLD, &status);
                // MPI_Get_count(&status, MPI_INT, &recv_count);
                // For now, assume fixed size
                recv_count = local_size;

                if (temp_buffer != NULL && recv_count > 0) {
                    MPI_Recv(temp_buffer, recv_count, MPI_INT, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // Place received data into the correct part of global_matrix
                    for (int i = 0; i < local_rows; i++) {
                        for (int j = 0; j < local_cols; j++) {
                            int global_i = r * local_rows + i;
                            int global_j = c * local_cols + j;
                            if (global_i < N && global_j < N) { // Boundary check
                                global_matrix[global_i * N + global_j] = temp_buffer[i * local_cols + j];
                            }
                        }
                    }
                } else if (recv_count > 0) { // temp_buffer is NULL but we expected data
                     fprintf(stderr, "Root: temp_buffer is NULL but expected to receive data from rank %d\n", source_rank);
                }
            }
        }
        if (temp_buffer) free(temp_buffer);

    } else {
        // Non-root processes send their local_data to root
        if (distrib->local_data != NULL && local_size > 0) {
            MPI_Send(distrib->local_data, local_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
    return 0;
}