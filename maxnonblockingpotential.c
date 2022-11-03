#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define Nx 10
#define varepsilon 0.001

//run to get rid of segfault: export OMPI_MCA_btl=self,tcp

void initialize(int argc, char *argv[] );

int rank, size;
int main(int argc, char *argv[]) {
    int Ny = Nx;    

    // MPI_Init(&argc, &argv);
    initialize(argc,argv);
    double t1 = MPI_Wtime();
    MPI_Status status;
    MPI_Request request;
    
    // CART TOPOL CREATED
    MPI_Comm comm1;
    int dimvec[2], periodvec[2], reorder;
    int id;
    int dimensions;
    dimensions = 2;
    dimvec[0] = 4; dimvec[1] = size/dimvec[0];
    periodvec[0] = 0; periodvec[1] = 0; 
    reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD,dimensions,dimvec,periodvec, reorder, &comm1);
    int coord[2];
    int findrank;
    int findcoord;
    int maxdims = 2;
    int up, down, right, left;
    for (findcoord =0; findcoord <size; findcoord++){
        MPI_Cart_coords(comm1, findcoord, maxdims, coord);
        MPI_Cart_rank(comm1, coord, &findrank);
        if (rank==findcoord){
            MPI_Cart_shift(comm1, 0, 1, &left, &right);
            MPI_Cart_shift(comm1, 1, 1, &up, &down);
            printf("%d:coord = (%3d,%3d); up=%3d; down=%3d; left=%3d; right=%3d \n",rank,coord[0],coord[1], up, down, left, right);
        }
    }

    double dx = 1.0/(Nx*dimvec[1]);
    double dy = 1.0/(Ny*dimvec[0]);

    // ALLOCATE MATRIX
    double(*locgrid)[Nx] = malloc(sizeof(double[Ny][Nx]));
    if(locgrid == NULL){
        printf("ERRRORR MALLOC\n");
    }
    for (int i =0; i<Ny; i++){
        for (int j = 0; j<Nx; j++){
            locgrid[i][j] = 1; 
        }
    }
    double(*nextlocgrid)[Nx] = malloc(sizeof(double[Ny][Nx]));
    if(nextlocgrid == NULL){
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<Nx-1; j++){
            nextlocgrid[i][j] = 0;
        }
    }

    // alloc boundary local buffers
    double *sendleft;
    sendleft = malloc((Ny/dimvec[1])*sizeof(double)); 
    double *recvleft;
    recvleft = malloc((Ny/dimvec[1])*sizeof(double));
    double *sendright;
    sendright = malloc((Ny/dimvec[1])*sizeof(double));
    double *recvright;
    recvright = malloc((Ny/dimvec[1])*sizeof(double));

    double *sendup;
    sendup = malloc((Nx/dimvec[0])*sizeof(double));
    double *recvup;
    recvup = malloc((Nx/dimvec[0])*sizeof(double));
    double *senddown;
    senddown = malloc((Nx/dimvec[0])*sizeof(double));
    double *recvdown;
    recvdown = malloc((Nx/dimvec[0])*sizeof(double));

    // Apply boundaries global buffer

    double *leftbdry;
    leftbdry = malloc(Ny*dimvec[1]*sizeof(double));
    double *rightbdry;
    rightbdry = malloc(Ny*dimvec[1]*sizeof(double));

    double *topbdry;
    topbdry = malloc(Nx*dimvec[0]*sizeof(double));
    double *botbdry;
    botbdry = malloc(Nx*dimvec[0]*sizeof(double));

    

    for (int i = 0; i < Nx*dimvec[0]-1; i++){
        botbdry[i] = 0;
        topbdry[i] = 0;
    }
    for (int i = 0; i < Ny*dimvec[1]-1; i++){
        leftbdry[i] = 0;
        rightbdry[i] = 0;
    }

    for (int i = 0.4*Nx; i < 0.6 * Nx * dimvec[0]; i++){
        botbdry[i] = 1;
    }



    int xcoord;
    int ycoord;
    // while (curr != 7){
    //     if (down == 0 || prev == left){
    //         for (int i = 0.6*Nx; i<Nx; i++){ // heat bot
    //                 locgrid[Ny-1][i] = 1;
    //         }
    //         prev = left;
    //     }
    // }
    // fill bdry

    xcoord = 0; // top bdry
    int xctr = 0;
    while (xctr < dimvec[1]){
        if (rank ==0){
            printf("rank %d: xcoord=%d\n",rank,xcoord);
        }
        xcoord = xcoord + dimvec[1];
        xctr++;
    }



    ycoord = 0; // left bdry
    int yctr = 0;
    while (yctr < dimvec[0]){
        if (rank ==0){
            printf("rank %d: ycoord=%d\n",rank,ycoord);
        }
        ycoord = ycoord +1;
        yctr++;
    }



    // for(ycoord = 0; ycoord<dimvec[1]; ycoord++){

    // }
    // for(xcoord = dimvec[0]-1; xcoord>-1; xcoord--){
        
    // }
    // for(ycoord = dimvec[1]-1; ycoord>-1; ycoord--){
        
    // }

    // if (rank ==0){
    //     printf("rank %d:\n",rank);
    //     for (int i =0; i<Ny; i++){
    //         for (int j = 0; j<Nx; j++){
    //             printf("%4.1f ", locgrid[i][j]);
    //         }
    //     printf("\n");
    //     }  
    // }




    int k = 0; 
    double *maxDelta;
    maxDelta = malloc(size*sizeof(double));
    maxDelta[0] = 1;
    int xstart=0;
    int ystart=0;
    double t2;
    int *breakflag;
    breakflag = malloc(2*sizeof(int));
    // printf("rank %d, finished malloc, starting loop\n",rank);
    while(1){
        k++;




        // find maxdelta local
        
        breakflag[0] = 0;
        if (rank == 0){
            if (k==100){
                breakflag[0]=1;
            }
            for (int i = 1; i<size; i++){
                // MPI_Recv(&maxDelta2,1,MPI_DOUBLE, i,1234,comm1,&status);
                // if (maxDelta[0]>maxDelta[1]){
                //     // if (maxDelta < varepsilon){
                //         breakflag[0] = 1;
                //         printf("%f maxdelta at rank: %d, breakflag = %d\n", maxDelta[0],rank,breakflag[0]);
                //     // }
                // }
                // if (maxDelta[1] > maxDelta[0]){
                //     maxDelta[0] = maxDelta[1];
                //     // if (maxDelta < varepsilon){
                //         breakflag[0] = 1;
                //         printf("%f maxdelta at rank: %d\n", maxDelta[0],rank);
                //     // }
                // }
            }
            // if (k%1000==0){
            //     // printf("%d iter, maxdelta = %f\n", k, maxDelta);
            // }
        }
    
        MPI_Bcast(breakflag, 1, MPI_INT, 0, comm1);
        for (int i = 1; i < Ny-1; i++){
            for (int j = 1; j < Nx-1; j++){
                locgrid[i][j] = nextlocgrid[i][j];
            }
        }
        // printf("%d , rank: %d\n", *breakflag, rank);
        if (breakflag[0] == 1){
            // if (rank == 0){
            //     t2 = MPI_Wtime();
            //     printf("my rank is: %d\nPoisson eq'n CONVERGED with Jacobi method\nunder varepsilon = %f, maxerror = %f \ntime taken = %12.8f seconds\n", 
            //     rank,varepsilon, maxDelta[0], t2-t1);
            // }
            break;
        }
    }
    // printf("%d: loopbroken\n",rank);
    MPI_Barrier(MPI_COMM_WORLD);
    free(leftbdry);
    free(rightbdry);
    free(topbdry);
    free(botbdry);
    free(breakflag);
    free(maxDelta);
    free(sendleft);
    free(recvleft);
    free(sendright);
    free(recvright);
    free(sendup);
    free(recvup);
    free(senddown);
    free(recvdown);
    free(locgrid);
    free(nextlocgrid);
    MPI_Comm_free(&comm1);
    // MPI_Request_free(&request);
    MPI_Finalize();
    return 0;
}

void initialize(int argc, char *argv[]) {
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
}