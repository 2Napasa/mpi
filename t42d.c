#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define N 400 // must be 400
#include <math.h>

int main(int argc, char *argv[]) {
    int Ny = N/2; 
    double t1 = MPI_Wtime();
    int rank, size;
    double dx = 1.0/(N*2);
    double dy = 1.0/(Ny*2);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    MPI_Request request, requestleft, requestright;
    


    // CREATING 2x2 cart TOPOL
    


    // CART TOPOL CREATED
    MPI_Comm comm1;
    int dimvec[2], periodvec[2], reorder;
    int id;
    int dimensions;
    // if (size !=6){
    //     printf("please run with 4 processess.\n"); fflush(stdout);
    //     MPI_Abort(MPI_COMM_WORLD,1);
    // }
    dimensions = 2;
    dimvec[0] = 2; dimvec[1] = 2;
    periodvec[0] = 0; periodvec[1] = 0; // tut nnado 0 postavit! net perioda
    reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD,dimensions,dimvec,periodvec, reorder, &comm1);
    // if (rank==0){
    //         printf("\n 2x2 comm1\n");
    // }
    int coord[2];
    int findrank;
    int findcoord;
    int maxdims = 2;
    int up, down, right, left;
    for (findcoord =0; findcoord <size; findcoord++){
        MPI_Cart_coords(comm1, findcoord, maxdims, coord);
        MPI_Cart_rank(comm1, coord, &findrank);
        // if (rank == 0){
        //     printf("for rank:%d found coord:(%d; %d); for coord found rank:%d\n", findcoord, coord[0],coord[1], findrank);
        // } 
        if (rank==findcoord){
            MPI_Cart_shift(comm1, 0, 1, &left, &right);
            MPI_Cart_shift(comm1, 1, 1, &up, &down);
            printf("%d: up=%2d; down=%2d; left=%2d; right=%2d\n", rank, up, down, left, right);
        }
    }


    // ALLOCATE MATRIX

    double(*locgrid)[N] = malloc(sizeof(double[Ny][N]));

    if(locgrid == NULL){
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<N-1; j++){
            locgrid[i][j] = 0; // IV fill 
        }
    }

    double(*nextlocgrid)[N] = malloc(sizeof(double[Ny][N]));
    if(nextlocgrid == NULL){
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

    double *sendup;
    sendup = malloc(N*size*sizeof(double));
    double *recvup;
    recvup = malloc(N*size*sizeof(double));
    double *senddown;
    senddown = malloc(N*size*sizeof(double));
    double *recvdown;
    recvdown = malloc(N*size*sizeof(double));

    double varepsilon = 0.000001;
    int k = 0; 
    double maxDelta;
    int xstart=0;
    double t2;
    while(1){
        k++;
        // BOUNDARIES DISABLED// // BOUNDARIES DISABLED// // BOUNDARIES DISABLED// // BOUNDARIES DISABLED// // BOUNDARIES DISABLED
        if (rank ==0){
            for (int i = 0*Ny; i<=.6*Ny; i++){ // heat 1 left
                locgrid[i][0] = 1;
            }
            
            for (int i = 0.6*N; i<N; i++){ // heat bot
                locgrid[0][i] = 1;
            }
        }
        if (rank == 2){
            for (int i = 0; i<0.6*N; i++){ // heat bot
                locgrid[0][i] = 1;
            }
            for (int i = 0; i<0.4*Ny; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 1 right
            }

            for (int i = 0.8*Ny; i<Ny; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 1 right
            }
            
            
        }
        if (rank == 1){
            for (int i = 0.2*Ny; i<=.6*Ny; i++){ // heat 2left 
                locgrid[i][0] = 1;
            }
        }
        if (rank == 3){
            for (int i = 0.0*Ny; i<0.2*Ny; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 2 right
            }
            for (int i = 0.4*Ny; i<0.8*Ny; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 3 right
            }
        }
        // BOUNDARIES DISABLED // // BOUNDARIES DISABLED// // BOUNDARIES DISABLED// // BOUNDARIES DISABLED// // BOUNDARIES DISABLED


        

        for (int i = 0; i < Ny; i++){
            for (int j = 0; j < N; j++){
                nextlocgrid[i][j]=locgrid[i][j];
            }
        }
        // solving local poisson // ENABLE POISSON // ENABLE LATER// ENABLE LATER// ENABLE LATER// ENABLE LATER// ENABLE LATER// ENABLE LATER

         // no need ystart since the f=x^2 given only depends on x

        if (rank == 2){
            xstart = N;
        }
        if (rank == 3){
            xstart = N;
        }
        for (int i = 1; i < Ny-1; i++){
            for (int j = 1; j < N-1; j++){
                // nextlocgrid[i][j] = 0.25 * (locgrid[i+1][j]+locgrid[i-1][j]+locgrid[i][j+1]+locgrid[i][j-1]);
                nextlocgrid[i][j] = ((locgrid[i+1][j]+locgrid[i-1][j])/dx/dx+(locgrid[i][j+1]+locgrid[i][j-1])/dy/dy+((j+xstart)*dx)*((j+xstart)*dx))/(2/dx/dx+2/dy/dy);
            }
        }
        // // COMMUNICATION TIME USING TOPOL
        
        // if (right != -2){
        //     for (int i = 0; i<N-2; i++){
        //         sendright[i] = locgrid[i+1][N-2];
        //     }
        // }
        // if (left!=-2){ 
        //     for (int i = 0; i<N-2; i++){
        //         sendleft[i] = locgrid[i+1][1];
        //     }
        // }
        // // printf("\n IM WORKING: %d\n", rank);
        // // }
        // // printf("%d: phantom arrays ready, starting grid exchange using topology, maxdelt = %9.6f\n",k, maxDelta);
        // if (left != -2){
        //     MPI_Isend(sendleft, N-1, MPI_DOUBLE,left,1234, dim1, &request);
        //     MPI_Irecv(recvleft, N-1, MPI_DOUBLE,left,1234,dim1, &request); 
        //     // MPI_Wait(&requestleft, &status);
        // }
        // // printf("\n IM WORKING: %d\n", rank);
        // if (right != -2){
        //     MPI_Isend(sendright, N-1, MPI_DOUBLE,right,1234,dim1, &request);
        //     // printf("\n %d: IM WORKING, my right neighbour is: %d\n",rank, right);
        //     MPI_Irecv(recvright, N-1, MPI_DOUBLE,right,1234,dim1, &request); 
        //     // MPI_Wait(&request, &status);
        // }
        // MPI_Wait(&request, &status);
        // // printf("%d: phantom grid exchange successful, maxdelt = %9.6f\n",k, maxDelta);
        
        // if (right!=-2){
        //     for (int i = 0; i<N-2; i++){
        //         locgrid[i+1][N-1] = recvright[i];
        //     }
        // }
        // if (left != -2){
        //     for (int i = 0; i<N-2; i++){
        //         locgrid[i+1][0] = recvleft[i];
        //     }
        // }
        // // COMMUNICATION ENDS

        // communication in X direction
        if (right != -2){
            for (int i = 0; i<Ny-2; i++){
                sendright[i] = locgrid[i+1][N-2];
            }
            MPI_Isend(sendright, Ny-1, MPI_DOUBLE,right,1234,comm1,&request);
            MPI_Irecv(recvright,Ny-1,MPI_DOUBLE, right,1234,comm1,&request);
            MPI_Wait(&request, &status);
            for (int i = 0; i<Ny-2; i++){
                locgrid[i+1][N-1] = recvright[i];
            }   
        }
        if (left != -2){
            for (int i = 0; i<Ny-2; i++){
                sendleft[i] = locgrid[i+1][1];
            }
            MPI_Isend(sendleft, Ny-1, MPI_DOUBLE,left,1234,comm1,&request);
            MPI_Irecv(recvleft,Ny-1,MPI_DOUBLE, left,1234,comm1,&request);
            MPI_Wait(&request, &status);
            for (int i = 0; i<Ny-2; i++){
                locgrid[i+1][0] = recvleft[i];
            }
            
        }

        // MPI_Barrier(MPI_COMM_WORLD);
        // COMMUNICATION IN X DIRECTION ENDS HERE


        // communication in Y direction STARTS
        if (down != -2){
            for (int i = 0; i<N-2; i++){
                senddown[i] = locgrid[Ny-2][i+1];
            }
            MPI_Isend(senddown, N-1, MPI_DOUBLE,down,1234,comm1,&request);
            MPI_Irecv(recvdown, N-1, MPI_DOUBLE,down,1234,comm1,&request);
            MPI_Wait(&request, &status);
            for (int i = 0; i<N-2; i++){
                locgrid[Ny-1][i+1] = recvdown[i];
            }   

            
        }
        if (up != -2){
            for (int i = 0; i<N-2; i++){
                sendup[i] = locgrid[1][i+1];
            }
            MPI_Isend(sendup, N-1, MPI_DOUBLE,up,1234,comm1,&request);
            MPI_Irecv(recvup, N-1, MPI_DOUBLE,up,1234,comm1,&request);
            MPI_Wait(&request, &status);
            for (int i = 0; i<N-2; i++){
                locgrid[0][i+1] = recvup[i];
            }
        }
        // proc 2-3
        // if (rank == 2){
        //     double *commbuffer;
        //     commbuffer = malloc(N*sizeof(double));
        //     for (int i = 0; i<N-2; i++){
        //         commbuffer[i] = locgrid[Ny-2][i+1];
        //     }
        //     MPI_Rsend(commbuffer, N-1, MPI_DOUBLE,3,1234,comm1);
        //     MPI_Recv(commbuffer,N-1,MPI_DOUBLE, 3,1234,comm1,&status);
        //     for (int i = 0; i<N-2; i++){
        //         locgrid[Ny-1][i+1] = commbuffer[i];
        //     }   
        // }
        // if (rank == 3){
        //     double *commbuffer;
        //     commbuffer = malloc(N*sizeof(double));
        //     MPI_Recv(commbuffer,N-1,MPI_DOUBLE, 2,1234,comm1,&status);
        //     for (int i = 0; i<N-2; i++){
        //         locgrid[0][i+1] = commbuffer[i];
        //     }
        //     for (int i = 0; i<N-2; i++){
        //         commbuffer[i] = locgrid[1][i+1];
        //     }
        //     MPI_Rsend(commbuffer, N-1, MPI_DOUBLE,2,1234,comm1);
        // }

        // COMMUNICATION IN Y DIRECTION ENDS HERE

        // CONVERGING LOCAL MAX
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
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,1234,comm1);
            // printf("%d: %d iter, maxdelta = %f\n", rank, k, maxDelta);
        }

        if (rank == 2){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,1234,comm1);
            // printf("%d: %d iter, maxdelta = %f\n", rank, k, maxDelta);
        }
        if (rank == 3){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,1234,comm1);
            // printf("%d: %d iter, maxdelta = %f\n", rank, k, maxDelta);
        }
        double maxDelta2, maxDelta3, maxDelta4;
        int breakflag=0;
        if (rank ==0 ){
            for (int i = 1; i<size; i++){
                MPI_Recv(&maxDelta2,1,MPI_DOUBLE, i,1234,comm1,&status);
                if (maxDelta>maxDelta2){
                    if (maxDelta < varepsilon){
                        breakflag = 1;
                        printf("%f \n", maxDelta);
                    }
                }
                if (maxDelta2 > maxDelta){
                    maxDelta = maxDelta2;
                    if (maxDelta < varepsilon){
                        breakflag = 1;
                        printf("%f \n", maxDelta);
                    }
                }
            }
            // if (k%1000==0){
            //     // printf("%d iter, maxdelta = %f\n", k, maxDelta);
            // }
        }
        MPI_Bcast(&breakflag, 1, MPI_INT, 0, comm1);
        for (int i = 1; i < Ny-1; i++){
            for (int j = 1; j < N-1; j++){
                locgrid[i][j] = nextlocgrid[i][j];
            }
        }
        // printf("%d , rank: %d\n",breakflag, rank);
        if (breakflag==1){
            if (rank == 0){
                t2 = MPI_Wtime();
                printf("my rank is: %d\nPoisson eq'n CONVERGED with Jacobi method\nunder varepsilon = %f, maxerror = %f \ntime taken = %12.8f seconds", 
                rank,varepsilon, maxDelta, t2-t1);
            }
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
            MPI_Finalize();
            return 0;
        }
    }

    free(locgrid);
    free(nextlocgrid);
    MPI_Finalize(); 
    return 0;
}