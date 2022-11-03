#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define N 200 // must be 200
#include <math.h>

int main(int argc, char *argv[]) {
    double t1 = MPI_Wtime();
    int Ny = N * 2;
    int rank, size;
    double dx = 1.0/(N*4);
    double dy = 1.0/(Ny);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(dim1, &size);
    MPI_Comm_rank(dim1, &rank);
    MPI_Status status;
    MPI_Request request, requestleft, requestright;

    

    // CREATE TOPOLOGY COMMUNICATOR
    if (rank == 0){
        printf("Creating 1x4 1D topology, shifts arer below:\n");
    }
    MPI_Comm dim1;
    int dimvec[2], periodvec[2], reorder;
    int id;
    int dimensions;
    dimensions = 1;
    dimvec[0] = 4; dimvec[1] = 1;
    periodvec[0] = 0; periodvec[1] = 0; // nonperiodic
    reorder = 1;
    MPI_Cart_create(dim1,dimensions,dimvec,periodvec, reorder, &dim1);
    int coord[2];
    int findrank;
    int findcoord;
    int maxdims = 2;
    int right, left;
    
    for (findcoord =0; findcoord <size; findcoord++){
        MPI_Cart_coords(dim1, findcoord, maxdims, coord);
        MPI_Cart_rank(dim1, coord, &findrank);
        if (rank==findcoord){
            MPI_Cart_shift(dim1, 0, 1, &left, &right);
            printf("%d: left=%3d; right=%3d\n", rank, left, right); // OUTCOMMENT // OUTCOMMENT // OUTCOMMENT // OUTCOMMENT // OUTCOMMENT // OUTCOMMENT // OUTCOMMENT // OUTCOMMENT // OUTCOMMENT 
        }
    }

    // COMMUNICATOR CREATED

    double(*locgrid)[N] = malloc (size * sizeof(double[Ny][N]));

    if(locgrid == NULL)
    {
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<N-1; j++){
            locgrid[i][j] = 0; // IV fill
        }
    }

    double(*nextlocgrid)[N] = malloc (size * sizeof(double[Ny][N]));
    if(nextlocgrid == NULL)
    {
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<N-1; j++){
            nextlocgrid[i][j] = 0; // IV fill
        }
    }

    double *sendleft;
    sendleft = malloc(N*size*sizeof(double));
    double *recvleft;
    recvleft = malloc(N*size*sizeof(double));
    double *sendright;
    sendright = malloc(N*size*sizeof(double));
    double *recvright;
    recvright = malloc(N*size*sizeof(double));

    double varepsilon = 0.0001;
    int k = 0; 
    double maxDelta = 1;
    double maxDelta2;
    int *breakflag;
    breakflag = malloc(size * sizeof(int));
    breakflag[0] = 0;
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
            for (int i = 0.6*N; i<N; i++){ // heat bot
                locgrid[0][i] = 1;
            } 
            // printf("rank%d: %d\n",rank,k);
            // for (int i =0; i<N; i++){
            //     for (int j = 0; j<N; j++){
            //         printf("%6.3f ", locgrid[i][j]);
            //     }
            //     printf("\n");
            // }
            
        }
        if (rank == 1){
            for (int i = 0.2*N; i<N; i++){ // heat bot
                locgrid[0][i] = 1;
            } 
            // printf("rank%d: %d\n",rank,k);
            // for (int i =0; i<N; i++){
            //     for (int j = 0; j<N; j++){
            //         printf("%6.3f ", locgrid[i][j]);
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
            // for (int i =0; i<Ny; i++){
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
        int xstart = N * rank;
        for (int i = 1; i < Ny-1; i++){
            for (int j = 1; j < N-1; j++){
                nextlocgrid[i][j] = ((locgrid[i+1][j]+locgrid[i-1][j])/dx/dx+(locgrid[i][j+1]+locgrid[i][j-1])/dy/dy+((xstart + j)*dx)*((xstart + j)*dx))/(2/dx/dx+2/dy/dy);
            }
        }

        // MPI_Barrier(dim1);
        // COMMUNICATION TIME USING TOPOL
        
        if (right != -2){
            for (int i = 0; i<N-2; i++){
                sendright[i] = locgrid[i+1][N-2];
            }
        }
        if (left!=-2){ 
            for (int i = 0; i<N-2; i++){
                sendleft[i] = locgrid[i+1][1];
            }
        }
        if (left != -2){
            MPI_Isend(sendleft, N-1, MPI_DOUBLE,left,1234, dim1, &request);
            MPI_Irecv(recvleft, N-1, MPI_DOUBLE,left,1234,dim1, &request); 
        }
        if (right != -2){
            MPI_Isend(sendright, N-1, MPI_DOUBLE,right,1234,dim1, &request);
            MPI_Irecv(recvright, N-1, MPI_DOUBLE,right,1234,dim1, &request); 
        }
        MPI_Wait(&request, &status);
        if (right!=-2){
            for (int i = 0; i<N-2; i++){
                locgrid[i+1][N-1] = recvright[i];
            }
        }
        if (left != -2){
            for (int i = 0; i<N-2; i++){
                locgrid[i+1][0] = recvleft[i];
            }
        }
        // COMMUNICATION ENDS
        
        
        // double *maxDelta;
        // maxDelta = malloc(size * sizeof(double));
        maxDelta = fabs(locgrid[0][0] - nextlocgrid[0][0]);
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
        
        // gathering local maxes
        if (rank == 1){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,555,dim1);
            // printf("%d: %d iter, maxdelta = %f\n", rank, k, maxDelta);
            
        }

        if (rank == 2){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,555,dim1);
            // printf("%d: %d iter, maxdelta = %f\n", rank, k, maxDelta);
        }
        if (rank == 3){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,555,dim1);
            // printf("%d: %d iter, maxdelta = %f\n", rank, k, maxDelta);
        }
        
        
        if (rank ==0 ){
            for (int i = 1; i<size; i++){
                MPI_Recv(&maxDelta2,1,MPI_DOUBLE, MPI_ANY_SOURCE,555,dim1,&status);
                if (maxDelta>maxDelta2){
                    if (maxDelta < varepsilon){
                        breakflag[0] = 1;
                        printf("%f terminating varepsilon \n", maxDelta);
                    }
                }
                if (maxDelta2 > maxDelta){
                    maxDelta = maxDelta2;
                    if (maxDelta < varepsilon){
                        breakflag[0] = 1;
                        printf("%f terminating varepsilon \n", maxDelta);
                    }
                }
            }
            
            
        }
        
        MPI_Bcast(breakflag, 1, MPI_INT, 0, dim1);
        // MPI_Barrier(dim1);
        for (int i = 1; i < Ny-1; i++){
            for (int j = 1; j < N-1; j++){
                locgrid[i][j] = nextlocgrid[i][j];
            }
        }
        // printf("%d , rank: %d\n",breakflag, rank);
        if (breakflag[0]==1){
            if (rank ==0 ){
                double t2 = MPI_Wtime();
                printf("my rank is: %d\nPoisson eq'n CONVERGED with Jacobi method\n maxerror = %f \ntime taken = %12.8f seconds\n", 
                    rank,varepsilon, t2-t1);
                printf("I am %d, finilazing my job here\n", rank);
            }
            // free(sendleft);
            // free(sendright);
            // free(recvleft);
            // free(recvright);
            // free(locgrid);
            // free(nextlocgrid);
            MPI_Finalize();
            return 0;
        }
        

    }
    // free(sendleft);
    // free(sendright);
    // free(recvleft);
    // free(recvright);
    // free(locgrid);
    // free(nextlocgrid);
    MPI_Finalize(); 
    return 0;
}