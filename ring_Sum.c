/***********************************************************************
* FILENAME :        ring_Sum.c
*
* DESCRIPTION :
*       Passing Messages from Proc to Proc. 
*
* NOTES :
*       Output to be handled by each process separately
* 	Sending Copies to rank 0 to be added to make output by single process
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
        int ProcNum, ProcRank, RecvRank;

        MPI_Status Status;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
        MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
        if (ProcNum < 1) {
                if (ProcRank == 0)  {
                        printf("use: \"mpicc -n 3 or more \"\n"); 
                }
                MPI_Finalize();
                return 0;
        }
        if (ProcRank == 0){
                printf("\nPerforming integer cycling test\n");
        }
        int sum = 0; 	// initial sum
        int recv[2];
        recv[1] = sum;
        for (int round = 1; round < 4; round++){
                if (ProcRank == 0){
                        int send[2];
                        send[0]=ProcRank;
                        send[1]=recv[1]+ProcRank;
                        int receiver = 1;
                        MPI_Rsend(&send, 2, MPI_INT,receiver,0,MPI_COMM_WORLD);
                        printf("Round%2d. ",round);
                        printf("I am%2d: ",ProcRank);
                        printf("sender:%2d; msg:%3d; to:%2d; status: sent\n",send[0],send[1],receiver);
                        MPI_Recv(&recv,2,MPI_INT,ProcNum-1,MPI_ANY_TAG,MPI_COMM_WORLD,&Status);
                        printf("Round%2d. ",round);
                        printf("I am%2d: ",ProcRank);
                        printf("sender:%2d; msg:%3d; to:%2d; status: received\n",recv[0],recv[1],ProcRank);
                }
                if (ProcRank != 0 & ProcRank != ProcNum-1){
                        MPI_Recv(&recv,2,MPI_INT,ProcRank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&Status);
                        printf("Round%2d. ",round);
                        printf("I am%2d: ",ProcRank);
                        printf("sender:%2d; msg:%3d; to:%2d; status: received\n",recv[0],recv[1],ProcRank);
                        int send[2];
                        send[0]=ProcRank;
                        send[1]=recv[1]+ProcRank;
                        int receiver = ProcRank+1;
                        MPI_Rsend(&send, 2, MPI_INT,receiver,0,MPI_COMM_WORLD);
                        printf("Round%2d. ",round);
                        printf("I am%2d: ",ProcRank);
                        printf("sender:%2d; msg:%3d; to:%2d; status: sent\n",send[0],send[1],receiver);
                }
                if (ProcRank == ProcNum-1){
                        MPI_Recv(&recv,2,MPI_INT,ProcRank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&Status);
                        printf("Round%2d. ",round);
                        printf("I am%2d: ",ProcRank);
                        printf("sender:%2d; msg:%3d; to:%2d; status: received\n",recv[0],recv[1],ProcRank);
                        int send[2];
                        send[0]=ProcRank;
                        send[1]=recv[1]+ProcRank;
                        int receiver = 0;
                        MPI_Rsend(&send, 2, MPI_INT,receiver,0,MPI_COMM_WORLD);
                        printf("Round%2d. ",round);
                        printf("I am%2d: ",ProcRank);
                        printf("sender:%2d; msg:%3d; to:%2d; status: sent\n",send[0],send[1],receiver);
                }
        }
        MPI_Finalize();
        return 0;
}