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


int main(int argc, char *argv[]) {
    int i, j, k, rank, size, tag = 12345, out = 0;
	int n = 200;
	
	int a[n][n];
	int b[n][n];
	int c[n][n];
    int j_gen, i_gen;
	// generating random matrix  with values 0-10
    srand(time(NULL));
    for (i_gen = 0; i_gen<n; i_gen++){
        for(j_gen = 0; j_gen<n; j_gen++){
            a[i_gen][j_gen] = rand() % 10; 
		}
	}
	srand(time(NULL)+1);	
	for (i_gen = 0; i_gen<n; i_gen++){
        for(j_gen = 0; j_gen<n; j_gen++){
            b[i_gen][j_gen] = rand() % 10;
		}
	}
	
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int a_local[n/size][n],b_local[n];
	int c_local[n/size][n];
	double t1 = MPI_Wtime();
        
    MPI_Scatter(a, n*n/size, MPI_INT, a_local, n*n/size, MPI_INT,0,MPI_COMM_WORLD);

	// for (i=0; i<size; i++){
	// 	if (i == rank){
	// 		printf("myrank %d: ------------------------------------\n",rank);
	// 		printf("Matrix A_local=\n");
	// 		for (i_gen = 0; i_gen < n/size; i_gen++){
	// 			for (j_gen = 0; j_gen < n; j_gen++){
	// 					printf("%3d\t", a_local[i_gen][j_gen]);
	// 			}
	// 			printf ("\n");
	// 		}
	// 	}
	// }

    MPI_Bcast(b, n*n, MPI_INT, 0, MPI_COMM_WORLD);
	int row, col;
	// SUBMATRIX MULTIPLICATION
	int l;
	for (i = 0; i < n/size; i++){ // by columns in out
		for (j = 0; j< n; j++){ // by rows in vector multiplication
			c_local[i][j]=0;
			for (l = 0; l < n; l++){ // by columns in vector multiplication
				out = out + a_local[i][l] * b[l][j];
			}
			c_local[i][j] = out;
			out = 0;
		}
	}
    MPI_Gather(&c_local, n*n/size, MPI_INT, c, n*n/size, MPI_INT, 0, MPI_COMM_WORLD);    
    
    if (rank == 0)  {
		
		double t2 = MPI_Wtime();
		double parallel_time = t2-t1;
		int c_0[n][n];
		for(i=0;i<n;i++){
                for(j=0;j<n;j++){
					c_0[i][j]=0;
					for(k=0;k<n;k++){
						c_0[i][j]+=a[i][k]*b[k][j];
				}
			}
        }
		double t3 = MPI_Wtime();
		double single_time = t3 - t2;
		// printf("myrank %d: ------------------------------------\n",rank);
		// printf("Matrix A=\n");
		// for (row = 0; row < n; row++) {
        //     for (col = 0; col < n; col++) {
        //             printf("%4d\t", a[row][col]);
        //     }
        //     printf ("\n");
    	// }
		// printf("Matrix B=\n");
		// for (row = 0; row < n; row++) {
        //     for (col = 0; col < n; col++) {
        //             printf("%4d\t", b[row][col]);
        //     }
        //     printf ("\n");
    	// }
		// printf("Matrix C = A x B =\n");
		// for (row = 0; row < n; row++) {
        //     for (col = 0; col < n; col++) {
        //             printf("%4d\t", c[row][col]);
        //     }
        //     printf ("\n");
    	// }
		printf("\ntime taken: \n1. parallel = %f\n2. oneproc = %f\n", parallel_time, single_time);
	}
	MPI_Finalize();
}