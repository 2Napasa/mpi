#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define N 100
#define Ny 200
#include <math.h>

int main(int argc, char *argv[]) {
    double t1 = MPI_Wtime();
    int rank, size;
    double dx = 1.0/(N*4);
    double dy = 1.0/(N);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    // store some values in x and y here

    double(*locgrid)[N] = malloc (sizeof(double[N][N]));

    if(locgrid == NULL)
    {
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<N-1; i++){
        for (int j = 1; j<N-1; j++){
            locgrid[i][j] = 0; // IV fill
        }
    }

    double(*nextlocgrid)[N] = malloc (sizeof(double[N][N]));
    if(nextlocgrid == NULL)
    {
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<N-1; i++){
        for (int j = 1; j<N-1; j++){
            nextlocgrid[i][j] = 0; // IV fill
        }
    }
    double varepsilon = 0.001;
    int k = 0; 
    while(1){
        k++;
        // boundaries

        // 
        if (rank ==0){
            for (int i = 0*N; i<=.3*N; i++){ // heat 1 left
                locgrid[i][0] = 1;
            }
            for (int i = 0.6*N; i<=.8*N; i++){ // heat 2left
                locgrid[i][0] = 1;
            }
            for (int i = 0.6*N; i<N; i++){ // heat bot
                locgrid[0][i] = 1;
            } 
            
        }
        if (rank == 1){
            for (int i = 0.2*N; i<N; i++){ // heat bot
                locgrid[0][i] = 1;
            } 
            // printf("rank%d: %d\n",rank,k);
            // for (int i =0; i<N; i++){
            //     for (int j = 0; j<N; j++){
            //         printf("%4.1f ", locgrid[i][j]);
            //     }
            //     printf("\n");
            // }  
        }

        if (rank == 2){
            for (int i = 0; i<0.4*N; i++){ // heat bot
                locgrid[0][i] = 1;
            }
        }
        if (rank == 3){
                    // for (int i = 0; i<0.2*N; i++){
                    //     locgrid[0][i] = 1;                  // heater bot
                    // }
            for (int i = 0; i<0.2*N; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 1 right
            }
            for (int i = 0.4*N; i<0.6*N; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 2 right
            }
            for (int i = 0.7*N; i<0.9*N; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 3 right
            }
        }
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                nextlocgrid[i][j]=locgrid[i][j];
            }
        }
        // solving local poisson
        for (int i = 1; i < N-1; i++){
            for (int j = 1; j < N-1; j++){
                // nextlocgrid[i][j] = 0.25 * (locgrid[i+1][j]+locgrid[i-1][j]+locgrid[i][j+1]+locgrid[i][j-1]);
                nextlocgrid[i][j] = ((locgrid[i+1][j]+locgrid[i-1][j])/dx/dx+(locgrid[i][j+1]+locgrid[i][j-1])/dy/dy+((j+1)*dx)*((j+1)*dx))/(2/dx/dx+2/dy/dy);
            }
        }
        // printf("time for communication\n");
        if (rank == 0){
            double *commbuffer;
            commbuffer = malloc(N*sizeof(double));
            for (int i = 0; i<N-2; i++){
                commbuffer[i] = locgrid[i+1][N-2];
            }
            // printf("my rank: %d, i m sending  to neighb", rank);
            MPI_Rsend(commbuffer, N-1, MPI_DOUBLE,1,1234,MPI_COMM_WORLD);
            
            MPI_Recv(commbuffer,N-1,MPI_DOUBLE, 1,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<N-2; i++){
                locgrid[i+1][N-1] = commbuffer[i];
            }   
        }
        if (rank == 1){
            double *commbuffer;
            commbuffer = malloc(N*sizeof(double));
            MPI_Recv(commbuffer,N-1,MPI_DOUBLE, 0,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<N-2; i++){
                locgrid[i+1][0] = commbuffer[i];
            }
            for (int i = 0; i<N-2; i++){
                commbuffer[i] = locgrid[i+1][1];
            }
            MPI_Rsend(commbuffer, N-1, MPI_DOUBLE,0,1234,MPI_COMM_WORLD);




            MPI_Recv(commbuffer,N-1,MPI_DOUBLE, 2,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<N-2; i++){
                locgrid[i+1][0] = commbuffer[i];
            }
            for (int i = 0; i<N-2; i++){
                commbuffer[i] = locgrid[i+1][1];
            }
            MPI_Rsend(commbuffer, N-1, MPI_DOUBLE,2,1234,MPI_COMM_WORLD);

        }
        if (rank == 2){
            double *commbuffer;
            commbuffer = malloc(N*sizeof(double));

            for (int i = 0; i<N-2; i++){
                commbuffer[i] = locgrid[i+1][1];
            }
            MPI_Rsend(commbuffer, N-1, MPI_DOUBLE,3,1234,MPI_COMM_WORLD);
            MPI_Recv(commbuffer,N-1,MPI_DOUBLE, 3,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<N-2; i++){
                locgrid[i+1][0] = commbuffer[i];
            }

            for (int i = 0; i<N-2; i++){
                commbuffer[i] = locgrid[i+1][1];
            }
            MPI_Rsend(commbuffer, N-1, MPI_DOUBLE,1,1234,MPI_COMM_WORLD);

            MPI_Recv(commbuffer,N-1,MPI_DOUBLE, 1,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<N-2; i++){
                locgrid[i+1][0] = commbuffer[i];
            }

            
            
        }

        if (rank == 3){
            double *commbuffer;
            commbuffer = malloc(N*sizeof(double));
            MPI_Recv(commbuffer,N-1,MPI_DOUBLE, 2,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<N-2; i++){
                locgrid[i+1][0] = commbuffer[i];
            }
            for (int i = 0; i<N-2; i++){
                commbuffer[i] = locgrid[i+1][1];
            }
            MPI_Rsend(commbuffer, N-1, MPI_DOUBLE,2,1234,MPI_COMM_WORLD);
        }


        // double *maxDelta;
        // maxDelta = malloc(size * sizeof(double));
        double maxDelta = fabs(locgrid[0][0] - nextlocgrid[0][0]);
        // if (rank == 1){
        //     printf("maxdelta changed to: %f\n ",nextlocgrid[0][0]);
        // }
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                if (fabs(locgrid[i][j] - nextlocgrid[i][j])>maxDelta){
                    maxDelta = fabs(locgrid[i][j] - nextlocgrid[i][j]);
                    // printf("maxdelta changed to: %f\n ",maxDelta);
                }
                // printf("maxdelta changed to: %f\n ",fabs(locgrid[i][j] - nextlocgrid[i][j]));
            }
        }
        // gathering local maxes
        if (rank == 1){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,1234,MPI_COMM_WORLD);
        }

        if (rank == 2){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,1234,MPI_COMM_WORLD);
        }
        if (rank == 3){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,1234,MPI_COMM_WORLD);
        }
        double maxDelta2, maxDelta3, maxDelta4;
        int breakflag=0;
        if (rank ==0 ){
            for (int i = 1; i<size; i++){
                MPI_Recv(&maxDelta2,1,MPI_DOUBLE, i,1234,MPI_COMM_WORLD,&status);
                if (maxDelta>maxDelta2){
                    if (maxDelta < varepsilon){
                        breakflag = 1;
                        // printf("%f \n", maxDelta);
                    }
                }
                if (maxDelta2 > maxDelta){
                    maxDelta = maxDelta2;
                    if (maxDelta < varepsilon){
                        breakflag = 1;
                        // printf("%f \n", maxDelta);
                    }
                }
            }
        }
        MPI_Bcast(&breakflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
        for (int i = 1; i < N-1; i++){
            for (int j = 1; j < N-1; j++){
                locgrid[i][j] = nextlocgrid[i][j];
            }
        }
        // printf("%d , rank: %d\n",breakflag, rank);
        if (breakflag==1){
            if (rank == 0){
                double t2 = MPI_Wtime();
                printf("my rank is: %d\n CONVERGED under varepsilon = %f, niter = %d\n time taken = %12.8f seconds", 
                rank,varepsilon,k, t2-t1);
            }
            free(locgrid);
            free(nextlocgrid);
            MPI_Finalize();
            return 0;
        }

    }

    free(locgrid);
    MPI_Finalize(); 
    return 0;
}