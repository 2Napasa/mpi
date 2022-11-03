/***********************************************************************
* FILENAME :        integral_mySum.c
*
* DESCRIPTION :
*       Calculating integral parallel and single-core. Comparing res.
*       
*
* NOTES :
*       MPI_Reduce and MPI_Op_create used
* 
* AUTHOR :    Rauan Kelesbekov        START DATE :    28 Nov 20
* CONTACTS :  rauan-k@hotmail.com
* 
* CHANGES : None
*
*
************************************************************************/
#include <mpi.h>
#include <stdio.h>
#include <math.h>

// functions descriptions are at the end of the code
long double trapezoidal(long double left, long double right, int count, long double dx);    
long double f(long double x);
void mySum (long double *, long double *, int *, MPI_Datatype *);


int main(int argc, char** argv) {
    long double a = 4;
    long double b = 5;
    int global_n = 100000;
    
    long double dx = (b-a)/global_n;
    MPI_Init(NULL, NULL);
    MPI_Status status;
    int procrank;
    int procnums;
    int tag = 0;
    long double t1 = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    MPI_Comm_size(MPI_COMM_WORLD, &procnums);
    MPI_Op myOp;
    
    long double local_n          = global_n/procnums;                                   // intervals on each proc
    long double local_left       = a + dx * procrank * local_n;                         // left side of each proc assigning
    long double local_right      = local_left + dx * local_n;                           // righ side of each proc assigning
    long double local_trapezoid  = trapezoidal(local_left, local_right, local_n, dx);   // calc inte on each proc
    long double out;

    MPI_Op_create((MPI_User_function *)mySum, 1, &myOp);
    MPI_Reduce(&local_trapezoid, &out, 1, MPI_LONG_DOUBLE, myOp, 0, MPI_COMM_WORLD);
    // send recv
    
    if (procrank == 0) { // all nonzero proc sending svoi int to rang0
        printf("____________________________________________________________________________________________\n");
        printf("\nI am %d: final ans = %13.10Lf;\n",procrank, out);
        double t2 = MPI_Wtime();
        double parallel_time= t2-t1;
        printf("I am %d: multi-proc using Reduce(myOp) computation is over, answ = %13.10Lf\n------- time taken = %13.10f sec \n",procrank,out,parallel_time);
        printf("_______________________________________________\n");
        double t_3              = MPI_Wtime();
        long double out_oneproc = trapezoidal(a, b, global_n, dx);
        double t_oneproc       = MPI_Wtime()-t_3;
        printf("I am %d: single-proc computation is over, answ = %13.10Lf\n------- time taken = %13.10f sec \n",procrank,out_oneproc,t_oneproc);
        printf("____________________________________________________________________________________________\n");
    }
    
    MPI_Finalize();
}


long double trapezoidal(long double left, long double right, int count , long double dx) {
   long double out, x; 
   int i;
   out = (f(left) + f(right))/2.0;
   for (i = 1; i <= count-1; i++) {
      x = left + i*dx;
      out += f(x);
   }
   out = out*dx;
   return out;
} 

long double f(long double x) {
   return 1/(x+1)/(x+2);
} 

void mySum(long double *toreduce, long double *reduced, int *len, MPI_Datatype *dtype)
{
    int i;
    // reduced = 0;
    for ( i=0; i<*len; i++ ) 
        reduced[0] += toreduce[0];
        printf("\n%Lf \n", reduced[i]);
}