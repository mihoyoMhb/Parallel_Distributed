#ifndef DATA_DISTRIBUTION_V2_H
#define DATA_DISTRIBUTION_V2_H

#include <mpi.h>

typedef struct {
    int world_size;      // Total number of processes
    int world_rank;      // Rank in MPI_COMM_WORLD
    
    int p_r;             // Number of processes in row dimension
    int p_c;             // Number of processes in column dimension
    
    int my_grid_row;     // This process\'s row position in the grid
    int my_grid_col;     // This process\'s column position in the grid
    
    int N;               // Global matrix size (N×N)
    int local_rows;      // Number of rows in local block (N/p_r)
    int local_cols;      // Number of columns in local block (N/p_c)
    
    MPI_Comm row_comm;   // Row communicator
    MPI_Comm col_comm;   // Column communicator
    
    int *local_data;     // Pointer to local data block (flattened)
} DataDistribution;

/**
 * Initialize the data distribution structure
 * @param N Global matrix size (N×N)
 * @param distrib Pointer to distribution structure to initialize
 * @return 0 on success, error code otherwise
 */
int init_distribution(int N, DataDistribution *distrib);

/**
 * Set up the processor grid and create communicators
 * @param distrib Initialized distribution structure
 * @return 0 on success, error code otherwise
 */
int setup_processor_grid(DataDistribution *distrib); // Will be modified for V2

/**
 * Allocate memory for local data block
 * @param distrib Initialized distribution structure
 * @return 0 on success, error code otherwise
 */
int allocate_local_block(DataDistribution *distrib);

/**
 * Distribute global matrix data from root to all processes
 * @param global_data Pointer to global matrix data (on root process only, can be NULL on other processes)
 * @param distrib Distribution structure with grid information
 * @return 0 on success, error code otherwise
 */
int distribute_matrix(int *global_data, DataDistribution *distrib);

/**
 * Reads matrix dimensions and data from a CSV file (root only),
 * then distributes the data to all processes.
 * N in distrib struct should be set (e.g. via bcast) before calling this.
 * @param filename Name of the CSV file
 * @param distrib Distribution structure
 * @param global_matrix_ptr_for_root_out Pointer to store the allocated global matrix data on root
 * @return 0 on success, error code otherwise
 */
int read_csv_and_distribute_data(const char* filename, DataDistribution *distrib, int** global_matrix_ptr_for_root_out);

/**
 * Gathers the distributed matrix data from all processes to the root process.
 * The root process must have allocated space for the global_matrix.
 * @param distrib Distribution structure
 * @param global_matrix Pointer to the buffer on root to store the gathered matrix (N x N)
 * @return 0 on success, error code otherwise
 */
int gather_full_matrix(DataDistribution *distrib, int *global_matrix);

/**
 * Prints the local data block for the current process.
 * @param distrib Distribution structure
 * @param title A title to print before the block
 */
void print_local_block(const DataDistribution *distrib, const char* title);

/**
 * Convert local indices to global indices
 * @param local_row Local row index
 * @param local_col Local column index
 * @param global_row Pointer to store global row index
 * @param global_col Pointer to store global column index
 * @param distrib Distribution structure
 */
void local_to_global_indices(int local_row, int local_col, 
                           int *global_row, int *global_col, 
                           const DataDistribution *distrib);

/**
 * Convert global indices to local indices
 * @param global_row Global row index
 * @param global_col Global column index
 * @param local_row Pointer to store local row index
 * @param local_col Pointer to store local column index
 * @param proc_row Pointer to store processor row
 * @param proc_col Pointer to store processor column
 * @param distrib Distribution structure
 * @return 1 if index belongs to this process, 0 otherwise
 */
int global_to_local_indices(int global_row, int global_col,
                          int *local_row, int *local_col,
                          int *proc_row, int *proc_col,
                          const DataDistribution *distrib);

/**
 * Clean up resources used by the distribution
 * @param distrib Distribution structure
 */
void cleanup_distribution(DataDistribution *distrib);

#endif /* DATA_DISTRIBUTION_V2_H */ 