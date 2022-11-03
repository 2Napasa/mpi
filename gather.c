#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include </usr/include/mpi/mpi.h>

int main ( int argc, char **argv ) {

    int rank;
    int size;

    int lvotes;
    int *gvotes;

    MPI_Init ( &argc, &argv );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &size );

    if (rank == 0)
        gvotes = malloc(size * sizeof(int) );

    /* everyone sets their first lvotes element */
    lvotes = rank+4;

    /* Gather to process 0 */
    MPI_Gather(&lvotes, 1, MPI_INT, /* send 1 int from lvotes.. */
                gvotes, 1, MPI_INT, /* gather 1 int each process into gvotes */
               0, MPI_COMM_WORLD); /* ... to root process 0 */


    printf("P%d: %d\n", rank, lvotes);
    if (rank == 0) {
        printf("P%d: Gathered ", rank);
        for (int i=0; i<size; i++)
            printf("%d ", gvotes[i]);
        printf("\n");
    }

    if (rank == 0)
        free(gvotes);

    MPI_Finalize();

    return 0;
}