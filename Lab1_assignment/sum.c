#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(int argc, char *argv[]){
    int n=0;
    // Check if the number of arguments is correct
    if(2 != argc){
        printf("Usage: sum number_of_summands\n");
        return 1;
    }
    int num_steps = atoi(argv[1]);
    MPI_Init(&argc, &argv);
    int myid, numprocs;
    // Get the number of processes and the rank of this process
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    //Using main process to initialize the array
    // double *a = NULL;
    // if(myid == 0){
    //     a = (double *)malloc(num_steps * sizeof(double));
    //     if(a == NULL){
    //         printf("Error allocating memory\n");
    //         MPI_Abort(MPI_COMM_WORLD, 1);
    //     }
    //     for(int i = 0; i < num_steps; i++){
    //         a[i] = (double)i;
    //     }
    // }

    /*Since the num_steps may not be able to be divised by size*/
    // The number of elements each process will handle
    int elelments_per_process = num_steps / numprocs;
    int remaining_elements = num_steps % numprocs;

    //Here, we use the scatterv function to send the data to all processes
    int *sendcounts = (int *)malloc(numprocs * sizeof(int));
    int *displacements = (int *)malloc(numprocs * sizeof(int));
    if(sendcounts == NULL || displacements == NULL){
        printf("Error allocating memory\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // Initialize the sendcounts and displacements arrays
    // Inspied by the link: https://rookiehpc.org/mpi/docs/mpi_scatterv/index.html
    for(int i = 0; i < numprocs; i++){
        if(i < remaining_elements){
            sendcounts[i] = elelments_per_process + 1;
        }
        else{
            sendcounts[i] = elelments_per_process;
        }
        if(i == 0){
            // The first process has no displacement
            displacements[i] = 0;
        }
        else{
            displacements[i] = displacements[i - 1] + sendcounts[i - 1];
        }
    }


    // // Allocate memory for the local array
    // double *my_array = (double *)malloc(sendcounts[myid] * sizeof(double));
    // int my_array_size = sendcounts[myid];
    // if(my_array == NULL){
    //     printf("Error allocating memory\n");
    //     MPI_Abort(MPI_COMM_WORLD, 1);
    // }
    // // Scatter the data to all processes
    // MPI_Scatterv(a, sendcounts, displacements, MPI_DOUBLE, my_array, my_array_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // 本地生成数据，避免分发


    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD); // 同步所有进程
    start_time = MPI_Wtime();


    int my_start = displacements[myid];
    int my_count = sendcounts[myid];
    double *my_array = malloc(my_count * sizeof(double));
    for (int i = 0; i < my_count; i++) {
        my_array[i] = my_start + i;
    }
    
    // 计算本地和
    double my_sum = 0.0;
    for (int i = 0; i < my_count; i++) {
           my_sum += my_array[i];
    }

    // Using tree reduction to calculate the total sum
    double total_sum = my_sum, temp;
    //Starting from the leaf nodes, see notes
    // https://github.com/HenryLiu0/MPI-Global-Summation/blob/master/%E6%A0%91%E5%BD%A2%E9%80%9A%E7%94%A8/%E6%A0%91%E5%BD%A2%E9%80%9A%E7%94%A8.c
    // int remaining_processes = numprocs, half, remainder;
    // while(remaining_processes != 1){
    //     half = remaining_processes / 2;
    //     remainder = remaining_processes % 2;
    //     if(myid < half){
    //         MPI_Recv(&temp, 1, MPI_DOUBLE, myid + half + remainder, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         total_sum += temp;
    //     }else if(myid >= half + remainder && myid < remaining_processes){
    //         MPI_Send(&my_sum, 1, MPI_DOUBLE, myid - half - remainder, 0, MPI_COMM_WORLD);
    //         break;
    //     }
    //     remaining_processes = half + remainder;
    // }
    int distance = 2;
    unsigned int max_distance = 1 << ((int)ceil(log2(numprocs)));
    while(distance <= max_distance){
        if(myid % distance == 0){
            // Measure the time taken for the reduction
            int source = myid + distance / 2;
            if(source < numprocs){
                MPI_Recv(&temp, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                total_sum += temp;
            }
        }else if(myid % distance == distance / 2){
            // Send the result
            MPI_Send(&total_sum, 1, MPI_DOUBLE, myid - distance / 2, 0, MPI_COMM_WORLD);
            break; 
        }
        distance *= 2;     
    }




    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    // 仅进程0输出时间
    if (myid == 0) {
        printf("Total sum: %f\n", total_sum);
        printf("Proc=%d | Problem size=%d | Time=%.6f sec\n", 
               numprocs, num_steps, end_time - start_time);
    }
    // Free allocated memory
    free(my_array);
    free(sendcounts);
    free(displacements);
    // Finalize MPI
    MPI_Finalize();
    return 0;
}


