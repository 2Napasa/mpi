/***********************************************************************
* FILENAME :        rankcheck.c
*
* DESCRIPTION :
*       Checking rank of a proc. 
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
int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // rank of process
    MPI_Comm_size(MPI_COMM_WORLD, &size); // nu of processes
    int count;
    for (int i = 0; i < size; i++) {
        if (i == rank){
            printf("I am %d of %d \n",rank, size);
            MPI_Barrier(MPI_COMM_WORLD);
            
        } 
        count++;  
    }
    if (rank == 0){
        printf("\n%d\n",count);
    }
    MPI_Finalize();
    return 0;
}