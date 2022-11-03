#include <stdio.h>
#include <stdlib.h>
#define Nx 10
#define Ny 10
int main(int argc, char *argv[]){
    
    double(*locgrid)[Nx] = malloc(sizeof(double[Ny][Nx]));
    if(locgrid == NULL){
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<Nx-1; j++){
            locgrid[i][j] = 0; 
        }
    }

    double(*nextlocgrid)[Nx] = malloc(sizeof(double[Ny][Nx]));
    if(nextlocgrid == NULL){
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<Nx-1; j++){
            nextlocgrid[i][j] = 1;
        }
    }
    locgrid = nextlocgrid;
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<Nx-1; j++){
            nextlocgrid[i][j] = 2;
        }
    }

    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<Nx-1; j++){
            printf("%4.2f ", locgrid[i][j]);
        }
        printf("\n");
    }
}