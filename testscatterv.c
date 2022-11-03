//MPI_Gatherv(&len,  100, MPI_DOUBLE, len_gathered, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define SIZE 8

int main(int argc, char *argv[])
{
    int rank, size;     
    int *sendcounts;    
    int *displs;        

    int rem = (SIZE)%size; 
    int sum = 0;                
    double rec_buf[100];          

    // the data to be distributed
    double data[SIZE] = {0,1,2,3,4,5,6,7};

    

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    sendcounts = malloc(sizeof(int)*size);
    displs = malloc(sizeof(int)*size);

    // calculate send counts and displacements
    for (int i = 0; i < size; i++) {
        sendcounts[i] = (SIZE)/size-1;
        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }

        displs[i] = sum;
        sum += sendcounts[i];
    }
    sendcounts[0]=1;    
    sendcounts[1]=2;
    sendcounts[2]=2;
    sendcounts[3]=3;
    displs[0] = 0;
    displs[1] = 1;
    displs[2] = 3;
    displs[3] = 5;

    // print calculated send counts and displacements for each process
    if (0 == rank) {
        for (int i = 0; i < size; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
    }

    // divide the data among processes as described by sendcounts and displs
    MPI_Scatterv(&data, sendcounts, displs, MPI_DOUBLE, &rec_buf, 100, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // print what each process received
    printf("%d: ", rank);
    for (int i = 0; i < sendcounts[rank]; i++) {
        printf("%lf\t", rec_buf[i]);
    }
    printf("\n");

    MPI_Finalize();

    free(sendcounts);
    free(displs);

    return 0;
}