/***********************************************************************
* FILENAME :        gatherv_scatterrv.c
*
* DESCRIPTION :
*       finding min/max length from set of vectors
*       
* NOTES :
*       Implemented 2 solutions: 
*           1. Calculation of min/max locally, sending results to root.
*           2. Calculation of lengths locally, sending array of len to root, 
*                   then finding min/max at the root.
*       Compared results.
*       Used Scatterv Gatherv for communication. 
* 
* AUTHOR :    Rauan Kelesbekov        START DATE :    20 Nov 20
* CONTACTS :  rauan-k@hotmail.com
* 
* CHANGES : None
*
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#define N 1000

int main(int argc, char *argv[])
{
    int *sendcounts;    
    int *displs;                    
    double rec_x[N]; 
    double rec_y[N];       
    double x[N];
    double y[N];

    srand(time(NULL));
    for (int i_gen = 0; i_gen<N; i_gen++){
            x[i_gen] = rand() % N+1; // remove "i_gen;//" for random array gen
	}
	srand(time(NULL)+1);	
	for (int i_gen = 0; i_gen<N; i_gen++){
            y[i_gen] = rand() % N+1; // remove "N-i_gen;//" for random array gen
	}
    int rank; 
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    sendcounts  = malloc(sizeof(int)*N);
    displs      = malloc(sizeof(int)*N);
    sendcounts[0] = .22*N;
    sendcounts[1] = .24*N;
    sendcounts[2] = .26*N;
    sendcounts[3] = .28*N;
    displs[0] = 0;
    for (int k=1; k< size; k++){
        displs[k] = displs[k-1] + sendcounts[k-1];
    }
    if (0 == rank) {
        for (int i = 0; i < size; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
    }

    MPI_Scatterv(&x, sendcounts, displs, MPI_DOUBLE, &rec_x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&y, sendcounts, displs, MPI_DOUBLE, &rec_y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int j = 0; j < size; j++){
        if (j == rank){
            printf("%d: \n", rank);
            for (int i = 0; i < sendcounts[rank]; i++) {
                printf("%3d. (%4.2lf, %4.2lf)\t ", i+1, rec_x[i], rec_y[i]);
            }
            printf("\n");
            
        }
        MPI_Barrier(MPI_COMM_WORLD); // to print in order
    }
    double min = pow(rec_x[0]*rec_x[0]+rec_y[0]*rec_y[0],0.5);
    double xmin = rec_x[0];
    double ymin = rec_y[0];
    double max = pow(rec_x[0]*rec_x[0]+rec_y[0]*rec_y[0],0.5);
    double xmax = rec_x[0];
    double ymax = rec_y[0];
    double len[sendcounts[rank]];
    for (int i=0; i<sendcounts[rank]; i++){
        if (pow(rec_x[i]*rec_x[i]+rec_y[i]*rec_y[i],0.5)<min){
            min = pow(rec_x[i]*rec_x[i]+rec_y[i]*rec_y[i],0.5);
            xmin = rec_x[i];
            ymin = rec_y[i];
        }
        if (pow(rec_x[i]*rec_x[i]+rec_y[i]*rec_y[i],0.5)>max){
            max = pow(rec_x[i]*rec_x[i]+rec_y[i]*rec_y[i],0.5);
            xmax = rec_x[i];
            ymax = rec_y[i];
        }
        len[i] = pow(rec_x[i]*rec_x[i]+rec_y[i]*rec_y[i],0.5);
    }
    // gathering calculated lengths
    double *len_gathered;
    len_gathered = malloc(sizeof(double)*N);
    MPI_Gatherv(&len,  sendcounts[rank], MPI_DOUBLE, len_gathered, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // gathering calculated min/max and indicies
    double min_gathered[size];
    double x_min_gathered[size];
    double y_min_gathered[size];
    double max_gathered[size];
    double x_max_gathered[size];
    double y_max_gathered[size];
    int gatherdispl[size];
    gatherdispl[0] = 0;
    gatherdispl[1] = 1;
    gatherdispl[2] = 2;
    gatherdispl[3] = 3;
    int rcounts[size];
    rcounts[0] = 1;
    rcounts[1] = 1;
    rcounts[2] = 1;
    rcounts[3] = 1;
    MPI_Gatherv(&min,  1, MPI_DOUBLE, min_gathered,     rcounts, gatherdispl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&max,  1, MPI_DOUBLE, max_gathered,     rcounts, gatherdispl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&xmin, 1, MPI_DOUBLE, x_min_gathered,   rcounts, gatherdispl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&ymin, 1, MPI_DOUBLE, y_min_gathered,   rcounts, gatherdispl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&xmax, 1, MPI_DOUBLE, x_max_gathered,   rcounts, gatherdispl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&ymax, 1, MPI_DOUBLE, y_max_gathered,   rcounts, gatherdispl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0){
        printf("\n_______________\n%d: gathered lengths (finding max/min when printing out):\n", rank);
        double minlen = len_gathered[0];
        double maxlen = len_gathered[0];
        double xminlen = x[0];
        double yminlen = y[0];
        double xmaxlen = x[0];
        double ymaxlen = y[0];
        for (int i=0; i<N; i++){
            if (len_gathered[i]<minlen){
                minlen = len_gathered[i];
                xminlen = x[i];
                yminlen = y[i];
            }
            if(len_gathered[i]>maxlen){
                maxlen = len_gathered[i];
                xmaxlen = x[i];
                ymaxlen = y[i];
            }
            printf("%6.2lf\t", len_gathered[i]);
        }
        printf("\n_______________\n%d: gathered minimums:\n", rank);
        for (int i = 0; i < 4; i++) {
            printf("%3d. %6.2lf at (%4.2lf, %4.2lf)\n", i+1, min_gathered[i],x_min_gathered[i],y_min_gathered[i]);
        }
        printf("\n");

        printf("___________\n%d: gathered maximums:\n", rank);
        for (int i = 0; i < 4; i++) {
            printf("%3d. %6.2lf at (%4.2lf, %4.2lf)\n", i+1, max_gathered[i],x_max_gathered[i],y_max_gathered[i]);
        }
        printf("\n");

        double min_global = min_gathered[0];
        double x_min_global = x_min_gathered[0];
        double y_min_global = y_min_gathered[0];
        double max_global = max_gathered[0];
        double x_max_global = x_max_gathered[0];
        double y_max_global = y_max_gathered[0];
        for (int i = 0; i < 4; i++) {
            if(min_global > min_gathered[i]){
                min_global = min_gathered[i];
                x_min_global = x_min_gathered[i];
                y_min_global = y_min_gathered[i];
            }
        }

        for (int i = 0; i < 4; i++) {
            if(max_global < max_gathered[i]){
                max_global = max_gathered[i];
                x_max_global = x_max_gathered[i];
                y_max_global = y_max_gathered[i];
            }
        }
        printf("local min/max to global - precalculating min/max on each proc separately\n");
        printf("min = %4.2lf; at (%4.2lf, %4.2lf)\n", min_global, x_min_global, y_min_global);
        printf("max = %4.2lf; at (%4.2lf, %4.2lf)\n", max_global, x_max_global, y_max_global);
        printf("len to min max - gathering lenghts to root finding min/max at root\n");
        printf("min = %4.2lf; at (%4.2lf, %4.2lf)\n", minlen, xminlen, yminlen);
        printf("max = %4.2lf; at (%4.2lf, %4.2lf)\n", maxlen, xmaxlen, ymaxlen);
        printf("if above results are the same, communication was completed successfully\n");
    }

    

    MPI_Finalize();

    free(sendcounts);
    free(displs);

    return 0;
}