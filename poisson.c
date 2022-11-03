#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#define n 1000000

int main ( int argc, char **argv ) {
    int rank;
    int size;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &size );
    MPI_Status status;
    double t1 = MPI_Wtime();
    double *global_grid_next;
    global_grid_next = malloc(size * n * sizeof(double));
    double *local_grid;
    double *local_grid_next;
    double *grid;
    grid = malloc(size * n * sizeof(double));
    double *global_grid;
    global_grid = malloc(size * n * sizeof(double));
    double left = 1; 
    double right = 0;
    grid[0] = left; // boundaries
    grid[n-1] = right; 
    // gen zero arrray
    for (int i_gen = 1; i_gen<n-1; i_gen++){
            grid[i_gen] = 0; 
            global_grid[i_gen] = grid[i_gen];
	}

    int i;
    double epsilon = 0.0001;
    double delta = 1;
    double gleft;
    double gright;
    int ctr = 0;
    while (1) {
        local_grid = malloc(size * (n/size+2) * sizeof(double));
        local_grid_next = malloc(size * (n/size+2) * sizeof(double));

        MPI_Scatter(grid, n/size, MPI_DOUBLE, local_grid, n/size, MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        // phantom array
        double *phantom_right;
        phantom_right = malloc(size * size * sizeof(double));
        for (int gr = 0; gr < size; gr++){
            phantom_right[gr] = grid[n/size * (gr+1)];
        }
        double *phantom_left;
        phantom_left = malloc(size * size * sizeof(double));
        phantom_left[0]=0;
        for (int gl = 1; gl < size; gl++){
            phantom_left[gl] = grid[n/size * gl-1];
            
        }
        //send phanntom array
        MPI_Scatter(phantom_left, 1, MPI_DOUBLE, &gleft, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Scatter(phantom_right, 1, MPI_DOUBLE, &gright, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        free(phantom_left);
        free(phantom_right);

        if (rank!=0){
            //shift elements to the right
            for (int k = n/size+1; k > 0; k--){        
                local_grid[k]=local_grid[k-1];
            }
            // phantom cell to the left 
            local_grid[0] = gleft;
            
        }
        if (rank!=size-1 && rank!=0){
            // phantom cell to the right
            local_grid[n/size+1] = gright;
            // printf("%d. phantom left received = %lf\n",rank, gleft);
            // printf("%d. phantom right received = %lf\n",rank, gright);
        }
        if (rank==0) {
            local_grid[n/size] = gright;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        // SOLVING LOCAL POISSON
        if (rank != 0 && rank != size -1) {
            for (int c = 1; c < n/size+1; c++){
                local_grid_next[c] = (local_grid[c-1] + local_grid[c+1])/2;
            }
        }
        else{
            for (int c = 1; c < n/size; c++){
                local_grid_next[c] = (local_grid[c-1] + local_grid[c+1])/2;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // reduce arrays to length n/size
        if (rank == 0){ 
            local_grid_next[0]=left;  // boundary left
        } 
        if (rank == size - 1){
            for (int clr = 0; clr<n/size+2; clr++){
                local_grid_next[clr]=local_grid_next[clr+1];
            }
            local_grid_next[n/size-1]=right; // boundary right
            // printf("\n right in last proc = %lf \n\n",local_grid_next[n/size-1]);
        }
        if(rank != size - 1 && rank !=0){
            for (int clr = 0; clr<n/size; clr++){
                local_grid_next[clr]=local_grid_next[clr+1]; 
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Gather(local_grid_next, n/size, MPI_DOUBLE, 
                grid, n/size, MPI_DOUBLE, 
                0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        
        //check for delta
        int breakflag = 0;
        if (rank == 0){
            double max_delta=0;
            for (int sv = 0; sv<n; sv++){
                if (fabs(global_grid[sv] - grid[sv]) > max_delta){
                    max_delta = fabs(global_grid[sv] - grid[sv]);   
                }
            }
            // printf("%d. max delta = %lf\n",ctr,max_delta);
            delta = max_delta;
            if (delta<epsilon){
                breakflag = 1;
            }
        }
        MPI_Bcast(&breakflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (breakflag==1){
            if (rank == 0){
                printf("my rank is: %d------------\n grid = ", rank);
                // for (int i = 0; i < n; i++) {
                //         printf("%5.2lf", grid[i]);
                // }
                printf ("\n");
                double t2 = MPI_Wtime();
                printf ("time = %lf\n",t2-t1);
            }
            free(local_grid);
            free(local_grid_next);
            free(global_grid_next);
            free(global_grid);
            free(grid);
            MPI_Finalize();
            return 0;
        }
        for (int sv = 0; sv<n; sv++){
            global_grid[sv] = grid[sv];
        }
        free(local_grid);
        free(local_grid_next);
        ctr++;
    }
    free(global_grid_next);
    free(global_grid);
    free(grid);
    MPI_Finalize();
    return 0;
}