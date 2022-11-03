#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int argc, char** argv) {
	MPI_Status status;
	MPI_Init(NULL, NULL);
	int ProcNum, ProcRank ;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum); 
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (ProcNum != 2) {
		fprintf(stderr, "ERRORERRORERROR");
		MPI_Abort(MPI_COMM_WORLD, 666 );
	}
	else {
		int itt = 0;
		int AddRank = (ProcRank +1 )%2; 
		int N = 100;
		while (itt < N) {
			if (ProcRank == itt % 2) { 
				int buffer = 1234; 
				itt++;
				MPI_Rsend( &itt, 1, MPI_INT, AddRank, 0, MPI_COMM_WORLD);
				printf("%d\t sent \t\t%d at \t%d'th iter\n",ProcRank, AddRank,itt);
 			}	 
			else {
				MPI_recv(&itt, 1, MPI_INT, AddRank, 0, MPI_COMM_WORLD, &status);
				printf("%d\t received \t%d at \t%d'th iter\n",ProcRank, AddRank, itt); 
 			}
		}
	}
 	MPI_Finalize();
}
