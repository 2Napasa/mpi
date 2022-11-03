/***********************************************************************
* FILENAME :        oddproc.c
*
* DESCRIPTION :
*       Working with odd proc's. 
*
* NOTES :
*       I/O to be handled by single process
* 
* AUTHOR :    Rauan Kelesbekov        START DATE :    7 Oct 20
* CONTACTS :  rauan-k@hotmail.com
* 
* CHANGES : None
*
*
************************************************************************/
#include <stdio.h>
#include <mpi.h>
#include <math.h>
int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // rank of process
    MPI_Comm_size(MPI_COMM_WORLD, &size); // nu of processes
    if (rank % 2 == 1)
        printf("I am %d of %d\n", rank, size);
    MPI_Finalize();
    return 0;
}