#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char*argv[]){
    int rank;
    MPI_Status status;

    struct{
        int x, y, z;
    } point;

    MPI_Datatype point_type;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Create a derived data type
    MPI_Type_contiguous(3, MPI_INT, &point_type);
    MPI_Type_commit(&point_type);
    // Initialize the data
    point.x = rank;
    point.y = rank + 1;
    point.z = rank + 2;
    
    printf("Rank %d: point = (%d, %d, %d)\n", rank, point.x, point.y, point.z);

    if(rank == 0){
        // Send the data to rank 1
        MPI_Send(&point, 1, point_type, 1, 0, MPI_COMM_WORLD);
        printf("Rank %d: Sent point = (%d, %d, %d)\n", rank, point.x, point.y, point.z);
    }else if(rank == 1){
        // Receive the data from rank 0
        MPI_Recv(&point, 1, point_type, 0, 0, MPI_COMM_WORLD, &status);
        printf("Rank %d: Received point = (%d, %d, %d)\n", rank, point.x, point.y, point.z);
    }
    // Free the derived data type
    MPI_Type_free(&point_type);
    MPI_Finalize();
    return 0;
}