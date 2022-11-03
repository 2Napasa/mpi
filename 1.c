#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define N 800
#define Ny 400
#include <math.h>

int main(int argc, char *argv[]) {
    double t1 = MPI_Wtime();
    int rank, size;
    double dx = 1.0/(N);
    double dy = 1.0/(Ny);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    // store some values in x and y here

    double(*locgrid)[N] = malloc (sizeof(double[Ny][N]));

    if(locgrid == NULL)
    {
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<N-1; j++){
            locgrid[i][j] = 0; // IV fill
        }
    }

    double(*nextlocgrid)[N] = malloc (sizeof(double[Ny][N]));
    if(nextlocgrid == NULL)
    {
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<N-1; j++){
            nextlocgrid[i][j] = 0; // IV fill
        }
    }
    double varepsilon = 0.000001;
    int k = 0; 
    
    while(1){
        k++;
        // boundaries
        if (rank ==0){
            for (int i = 0*Ny; i<=.3*Ny; i++){ // heat 1 left
                locgrid[i][0] = 1;
            }
            for (int i = 0.6*Ny; i<=.8*Ny; i++){ // heat 2left
                locgrid[i][0] = 1;
            }
            for (int i = 0.3*N; i<0.6*N; i++){ // heat bot
                locgrid[0][i] = 1;
            }

            for (int i = 0; i<0.2*Ny; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 1 right
            }
            for (int i = 0.4*Ny; i<0.6*Ny; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 2 right
            }
            for (int i = 0.7*Ny; i<0.9*Ny; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 3 right
            }
            // printf("rank%d: %d\n",rank,k);
            // for (int i =0; i<N; i++){
            //     for (int j = 0; j<N; j++){
            //         printf("%4.1f ", locgrid[i][j]);
            //     }
            //     printf("\n");
            // }  
        }
        for (int i = 0; i < Ny; i++){
            for (int j = 0; j < N; j++){
                nextlocgrid[i][j]=locgrid[i][j];
            }
        }
        // solving local poisson
        for (int i = 1; i < Ny-1; i++){
            for (int j = 1; j < N-1; j++){
                // nextlocgrid[i][j] = 0.25 * (locgrid[i+1][j]+locgrid[i-1][j]+locgrid[i][j+1]+locgrid[i][j-1]);
                nextlocgrid[i][j] = ((locgrid[i+1][j]+locgrid[i-1][j])/dx/dx+(locgrid[i][j+1]+locgrid[i][j-1])/dy/dy+((j+1)*dx)*((j+1)*dx))/(2/dx/dx+2/dy/dy);
            }
        }

        // if (rank == 0){
        //     double *commbuffer;
        //     commbuffer = malloc(N*sizeof(double));
        //     for (int i = 0; i<N-2; i++){
        //         commbuffer[i] = locgrid[i+1][N-2];
        //     }
        //     MPI_Rsend(commbuffer, N-1, MPI_DOUBLE,1,1234,MPI_COMM_WORLD);
        //     MPI_Recv(commbuffer,N-1,MPI_DOUBLE, 1,1234,MPI_COMM_WORLD,&status);
        //     for (int i = 0; i<N-2; i++){
        //         locgrid[i+1][N-1] = commbuffer[i];
        //     }   
        // }
        // if (rank == 1){
        //     double *commbuffer;
        //     commbuffer = malloc(N*sizeof(double));
        //     MPI_Recv(commbuffer,N-1,MPI_DOUBLE, 0,1234,MPI_COMM_WORLD,&status);
        //     for (int i = 0; i<N-2; i++){
        //         locgrid[i+1][0] = commbuffer[i];
        //     }
        //     for (int i = 0; i<N-2; i++){
        //         commbuffer[i] = locgrid[i+1][1];
        //     }
        //     MPI_Rsend(commbuffer, N-1, MPI_DOUBLE,0,1234,MPI_COMM_WORLD);
        // }


        // double *maxDelta;
        // maxDelta = malloc(size * sizeof(double));
        double maxDelta = fabs(locgrid[0][0] - nextlocgrid[0][0]);
        // if (rank == 1){
        //     printf("maxdelta changed to: %f\n ",nextlocgrid[0][0]);
        // }
        for (int i = 0; i < Ny; i++){
            for (int j = 0; j < N; j++){
                if (fabs(locgrid[i][j] - nextlocgrid[i][j])>maxDelta){
                    maxDelta = fabs(locgrid[i][j] - nextlocgrid[i][j]);
                    // printf("maxdelta changed to: %f\n ",maxDelta);
                }
                // printf("maxdelta changed to: %f\n ",fabs(locgrid[i][j] - nextlocgrid[i][j]));
            }
        }
        // printf("maxdelta curr at %d iter: %f\n ",k,maxDelta);
        // gathering local maxes
        // if (rank == 1){
        //     MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,1234,MPI_COMM_WORLD);
        // }
        
        int breakflag=0;
        if (maxDelta<varepsilon){
            breakflag=1;
        }
        for (int i = 1; i < Ny-1; i++){
            for (int j = 1; j < N-1; j++){
                locgrid[i][j] = nextlocgrid[i][j];
            }
        }
        
        // printf("%d , rank: %d\n",breakflag, rank);
        if (breakflag==1){
            if (rank == 0){
                double t2 = MPI_Wtime();
                printf("my rank is: %d\nPoisson eq'n CONVERGED with Jacobi method\nunder varepsilon = %f, maxerror = %f \ntime taken = %12.8f seconds", 
                rank,varepsilon, maxDelta, t2-t1);
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