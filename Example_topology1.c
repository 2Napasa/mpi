#include <mpi.h>
#include <stdio.h>
int main(int argc, char *argv[]){
    int rank, size;
    MPI_Comm ourLove;
    int dimvec[2] = {2,2};
    int periodvec[2] = {0,0};
    int reorder = 1;
    int dimensions = 2; 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Cart_create(MPI_COMM_WORLD,dimensions,dimvec,periodvec, reorder, &ourLove);
    int coord[2];
    int findrank;
    int findcoord;
    int maxdims = 2;
    for (findcoord =0; findcoord <size; findcoord++){
        MPI_Cart_coords(ourLove, findcoord, maxdims, coord);
        MPI_Cart_rank(ourLove, coord, &findrank);
        if (rank==findcoord){
            int up, down, right, left;
            MPI_Cart_shift(ourLove, 0, 1, &left, &right);
            MPI_Cart_shift(ourLove, 1, 1, &up, &down);
        }
    }
    MPI_Finalize();
} 



// MPI_Comm comm2;
//     dimensions = 2;
//     dimvec[0] = 4; dimvec[1] = 1;
//     periodvec[0] = 0; periodvec[1] = 0; // tut nnado 0 postavit! net perioda
//     reorder = 1;
//     MPI_Cart_create(MPI_COMM_WORLD,dimensions,dimvec,periodvec, reorder, &comm2);
//     maxdims = 2;

//     for (findcoord =0; findcoord <size; findcoord++){
//         MPI_Cart_coords(comm2, findcoord, maxdims, coord);
//         MPI_Cart_rank(comm2, coord, &findrank);
//         // if (rank == 0){
//         //     printf("for rank:%d found coord:(%d; %d); for coord found rank:%d\n", findcoord, coord[0],coord[1], findrank);
//         // } 
//         // if (rank==findcoord){
//         //     int up, down, right, left;
//         //     MPI_Cart_shift(ourLove, 0, 1, &left, &right);
//         //     MPI_Cart_shift(ourLove, 1, 1, &up, &down);
//         //     printf("%d: up=%d; down=%d; left=%d; right=%d\n", rank, up, down, left, right);
//         // }
//     }

