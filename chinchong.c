/***********************************************************************
* FILENAME :        pingpong.c
*
* DESCRIPTION :
*       Passing Messages from Proc to Proc. 
*
* NOTES :
*       I/O to be handled by single process
* 		Was willing to make a pp code for even nu of proc, then realized i only have 2 cores :-(
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
	if (ProcNum != 2) {
		if (ProcRank == 0)  {
			printf("use: \"mpicc -n 2 \"\n"); 
		}
		MPI_Finalize();
		return 0;
	}
	
    for (int j = 1; j<101; j++){
		if (ProcRank == 0){
				/* If you require that your messages are printed in a specific order, 
				you must send them all to one rank which can print them in whatever order you like. 
				As long as one rank does all of the printing, all of the messages will be ordered correctly.
				https://stackoverflow.com/questions/17570996/mpi-printing-in-an-order	*/	
				MPI_Rsend(&ProcRank, 1, MPI_INT,-2,0,MPI_COMM_WORLD);
				printf("|iter%3d | sent | %2d ||\n",j, ProcRank);
				MPI_Recv(&RecvRank, 1, MPI_INT,-2,MPI_ANY_TAG,MPI_COMM_WORLD,&Status);
				printf("|iter%3d |  got | %2d ||\n",j, RecvRank);	
		}
		else {
				MPI_Recv(&RecvRank, 1, MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&Status);
				MPI_Send(&ProcRank, 1, MPI_INT,0,0,MPI_COMM_WORLD);
		}
	}
	if (ProcRank == 0){
		printf("\nPerforming integer bouncing test\n");
	}
	int ball 	= 0; 	// bounces
	int hits0 	= 0; 	// hits by proc0
	int hits1 	= 0; 	// hits by proc1
	int hits1_ball[2]; 	// array to send hits and bounces
	for (int j = 1; j < 101; j++){
		if (ProcRank == 0){
				ball++;
				hits0++;
				// sending total nu of hits
				MPI_Send(&ball, 1, MPI_INT,-2,0,MPI_COMM_WORLD);
				printf("|iter%3d |  %d ping | ball bounces = %3d | Proc_0 hits %3d times |\n",j,ProcRank, ball,hits0);
				// receiving nu of hits by proc1 and a total nu of ball bounces
				MPI_Recv(&hits1_ball, 2, MPI_INT,-2,MPI_ANY_TAG,MPI_COMM_WORLD,&Status); 
				ball = hits1_ball[0]; // storing received number of hits
				printf("|iter%3d |  %d pong | ball bounces = %3d | Proc_1 hits %3d times |\n",j,ProcRank+1, hits1_ball[0], hits1_ball[1]);	
		}
		else {
				// receiving total nu of hits
				
				hits1++; ball++; // increasing hits on proc1 and a total nu of ball bounces
				hits1_ball[0] = ball; // preparing an array to send 
				hits1_ball[1] = hits1;
				// sending total nu of bounces and nu of hits by proc1
				MPI_Send(&hits1_ball, 2, MPI_INT,0,0,MPI_COMM_WORLD);
				MPI_Recv(&ball, 1, MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&Status);
		}
	}
    MPI_Finalize();
    return 0;
}