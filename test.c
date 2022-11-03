
/***********************************************************************
* FILENAME :        matMult_mpi.c
*
* DESCRIPTION :
*       Multiplicating matrix parallelly A x B = C
*       
*
* NOTES :
*       used Scatter Gather for communication. 
* 		
* 
* AUTHOR :    Rauan Kelesbekov        START DATE :    2 Nov 20
* CONTACTS :  rauan-k@hotmail.com
* 
* CHANGES : None
*
*
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#define n 4

int main(int argc, char *argv[]) {
    int i, j, k, rank, size, tag = 12345, out = 0;
	
	int a[n][n];
	int b[n][n];
    int iii, ooo;
    // srand(time(NULL));
    for(ooo = 0; ooo<n; ooo++)
        for(iii = 0; iii<n; iii++)
            a[ooo][iii] = iii+ooo+iii; //rand() % 5; [] //0     1     2     3
													//1     2     3     4
													//2     3     4     5
													//3     4     5     6
	// srand(time(NULL));	
	
	for(ooo = 0; ooo<n; ooo++)
        for(iii = 0; iii<n; iii++)
            b[ooo][iii] = iii+ooo+2*iii; //rand() % 5;
	
    int c[n][n];
	// if (rank == 0) {
	// 	printout("A = ", a);
	// 	printout("B = ", b);
	// }
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int a_local[n/size][n],b_local[n];
	int c_local[n][n/size];
	// if (rank == 0 ){
	double t1 = MPI_Wtime();
	// }
        
    // MPI_Scatter(a, n*n/size, MPI_INT, a_local, n*n/size, MPI_INT,0,MPI_COMM_WORLD);
	int row, col;
	int m = n/size;
	if (rank ==0){
		printf("myrank %d: ------------------------------------\n",rank);
		printf("Matrix A=\n");
		for (row = 0; row < n; row++) {
            for (col = 0; col < n; col++) {
                    printf("%3d\t", a[row][col]);
            }
            printf ("\n");
    	}
		printf("Matrix B=\n");
		for (row = 0; row < n; row++) {
            for (col = 0; col < n; col++) {
                    printf("%3d\t", b[row][col]);
            }
            printf ("\n");
    	}
		printf("\n");
		// printf("Matrix A_local=\n");
    	// for (row = 0; row < m; row++) {
        //     for (col = 0; col < n; col++) {
        //             printf("%3d\t", a_local[row][col]);
        //     }
        //     printf ("\n");
    	// }
	}
	// else {
	// 	printf("myrank %d: ------------------------------------ \n",rank);
	// 	printf("Matrix A_local=\n");
    // 	for (row = 0; row < m; row++) {
    //         for (col = 0; col < n; col++) {
    //                 printf("%3d\t", a_local[row][col]);
    //         }
    //     	printf ("\n");
    // 	}
	// }
	printf("\n");



	// MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Bcast(b, n*n, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);


    // MPI_Gather(c_local, n*n/size, MPI_INT, c, n*n/size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Finalize();
}