#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void mySum(int *toreduce, int *reduced, int *len, MPI_Datatype *dtype)
{
    int i;
    // reduced = 0;
    for ( i=0; i<*len; i++ ) 
        reduced[i] = reduced[i] + toreduce[i];
        printf("\n, reduced[%d] = %d; toreduce[%d] = %d \n", i, reduced[i], i, toreduce[i]);
}
int main(int argc, char **argv) {
    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *globaldata=NULL;
    int localdata;

    if (rank == 0) {
        globaldata = malloc(size * sizeof(int) );
        for (int i=0; i<size; i++)
            globaldata[i] = 2*i+1;

        printf("Processor %d has data: ", rank);
        for (int i=0; i<size; i++)
            printf("%d ", globaldata[i]);
        printf("\n");
    }

    MPI_Scatter(globaldata, 1, MPI_INT, 
                &localdata, 1, MPI_INT, 
                0, MPI_COMM_WORLD);

    printf("Processor %d has data %d\n", rank, localdata);
    localdata *= 2;
    printf("Processor %d doubling the data, now has %d\n", rank, localdata);

    // MPI_Gather(&localdata, 1, MPI_INT, 
    //             globaldata, 1, MPI_INT, 
    //             0, MPI_COMM_WORLD);

    MPI_Op myOp;
    MPI_Op_create((MPI_User_function *)mySum, 1, &myOp);
    MPI_Reduce(&localdata, globaldata, 1, MPI_INT, myOp, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Processor %d has data: ", rank);
        for (int i=0; i<size; i++)
            printf("%d ", globaldata[i]);
        printf("\n");
    }

    if (rank == 0)
        free(globaldata);

    MPI_Finalize();
    return 0;
}