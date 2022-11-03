/***********************************************************************
* FILENAME :        ring.c
*
* DESCRIPTION :
*       Non-Blocking passing  to the neighbour, getting from the neighbour. 
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
    int my_rank, ncpus;
    int left_neighbor, right_neighbor;
    int data_received = -1;
    int tag = 101;
    MPI_Status statSend, statRecv;
    MPI_Request reqSend, reqRecv;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
    left_neighbor = (my_rank-1+ncpus)%ncpus;
    right_neighbor = (my_rank+1)%ncpus;

    MPI_Isend(&my_rank,         1,MPI_INT, 
            left_neighbor,   tag, MPI_COMM_WORLD, &reqSend);
    MPI_Irecv(&data_received,   1,MPI_INT, 
            right_neighbor,  tag, MPI_COMM_WORLD, &reqRecv);
    MPI_Wait(&reqSend, &statSend);
    MPI_Wait(&reqRecv, &statRecv);
    printf("Among %d processes, process %d received from right neighbor: %d\n",
            ncpus, my_rank,data_received);
    MPI_Finalize();    
    return 0;

}