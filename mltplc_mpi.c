/* 
   A simple matrix-matrix multiplication example using 
   MPI_SCATTER, MPI_GATHER, MPI_BCAST.

   The example arrays are dynamically allocated arrays, and the
   the number of rows in the product matrix must be divisible by
   the number of processors.

   5 April 1998
   Revised September 26, 1998 
   Andrew C. Pineda
   Albuquerque High Performance Computing Center
   1601 Central Ave NE
   Albuquerque, NM 87131
   acpineda@arc.unm.edu
   $Id: matmat.c,v 1.1 1998/04/07 12:17:12 pineda Exp acpineda $

  */

/* 
   For programmers with no experience with C language programming.
   The C language has little in the way of built-in operators and
   functions. These instead are accessed via library calls. STDIO.H
   provides definitions for basic I/O functions, e.g. printf. STDLIB.H
   provides definitions for the functions malloc and free which are
   the functions used to allocate and deallocate memory from the
   system. MPI.H provides the definitions for MPI functions and
   constants including timing calls. 
*/ 

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h> 

/* 
   Allocate a pointer to a 2d array so that the data is stored in a
   contiguous piece of memory, so that MPI can distribute it
   correctly. 

   Quick C lesson for non-C programmers: Anything with the unary
   operator, * , in front of it is a pointer quantity. Something with
   two *'s in front of it is a pointer to a pointer, etc. The key to
   decoding what a pointer variable is to look at its declaration.
   For example, if the declaration of a variable b is
     int *b;
   Then, when *b appears in an expression, it is an int, and when b
   appears in an expression it is a pointer to an it. Pointer and
   array syntax are frequently used interchangably in C. In the example 
   above, if b points to a location in memory large enough to hold 10
   ints. Then,
   
   b[0] is the same as *b
   b[1] is the same as *(b+1)
   b[2] is the same as *(b+2)
   etc.

*/

double **allocate_array(int row_dim, int col_dim) 
{
  double **result;
  int i;

  /* Allocate an array of pointers to hold pointers to the rows of the
	 array */
  result=(double **)malloc(row_dim*sizeof(double *));

  /* The first pointer is in fact a pointer to the entire array */
  result[0]=(double *)malloc(row_dim*col_dim*sizeof(double));

  /* The remaining pointers are just pointers into this array, offset
	 in units of col_dim */
  for(i=1; i<row_dim; i++)
	result[i]=result[i-1]+col_dim;

  return result;
}

/* Deallocate a pointer to a 2d array */

void deallocate_array(double **array, int row_dim) 
{
  int i;
  /* Make sure all the pointers into the array are not pointing to
	 random locations in memory */
  for(i=1; i<row_dim; i++)
	array[i]=NULL;
  /* De-allocate the array */
  free(array[0]);
  /* De-allocate the array of pointers */
  free(array);
}

/* For this cooked up example, the result of the product is known in
   closed form. EXACT computes the (i,j)th element of the product of
   A and B from this expression. The calculation is done in double
   precision because the (i,j)th term may be larger than what an int
   can hold. */

double exact(int i, int j, int k) 
{
  double result;
  double ri, rj, rk;
  ri=(double) i;
  rj=(double) j;
  rk=(double) k;
  result=(ri*rj*rk)+((ri+rj)*rk*(rk-1.0)/2.0)+((rk-1.0)*rk*(2.0*rk-1.0)/6.0);
  return result;
}

/* The main body of the program. Note to non-C Programmers, every C
   program must have a function called MAIN. MAIN can receive command
   line arguments via ARGC and ARGV. */

int main( int argc, char **argv )
{
  /* Declarations of variables from serial version */
  double **a, **b, **c;
  int nrow_a, ncol_a, nrow_b, ncol_b, nrow_c, ncol_c;
  int i,j,k;

  /* Timing */
  double start_time, end_init, end_calc, end_time;

  /* Additional MPI variables */
  double **apart,  **cpart;
  int rank, size, ierr;
  int root=0;

  /* Execution starts here */

  /* Start up MPI communications */
  ierr=MPI_Init(&argc, &argv);
  ierr=MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  ierr=MPI_Comm_size(MPI_COMM_WORLD,&size);

  /* Only do I/O and timing on the root processor */
  if(rank==root) {

    printf("Enter rows of A, cols of A=rows of B, cols of B:\n");
    scanf("%d %d %d",&nrow_a, &ncol_a, &ncol_b);

    start_time=MPI_Wtime();
  }

  /* Wait for the root processor to do the I/0. */
  /* What would happen if the Barrier call is removed? */
  ierr=MPI_Barrier(MPI_COMM_WORLD);

  /* Broadcast the column and row dimensions. Note:
     This is not the most efficient way to do this.
     I.e., we really should package the 3 ints into
     an array and broadcast the array, but since this
     is a one time initialization we don't bother. */

  ierr=MPI_Bcast(&nrow_a,1,MPI_INT,root,MPI_COMM_WORLD);
  ierr=MPI_Bcast(&ncol_a,1,MPI_INT,root,MPI_COMM_WORLD);
  ierr=MPI_Bcast(&ncol_b,1,MPI_INT,root,MPI_COMM_WORLD);

  /* Now that everyone knows the array dimensions,
     allocate arrays */

  nrow_b=ncol_a;
  nrow_c=nrow_a;
  ncol_c=ncol_b;

  b=allocate_array(nrow_b,ncol_b); /* Everyone has a copy of this in 
                                      our simple MPI example */

  /* Only root needs to have all of product and right hand factor */
  if(rank==root) {
    a=allocate_array(nrow_a,ncol_a);
    c=allocate_array(nrow_c,ncol_c);
  }

  /* Everyone allocates space for their sub-blocks of A and C. */
  apart=allocate_array(nrow_a/size,ncol_a);
  cpart=allocate_array(nrow_c/size,ncol_c);

  /* Initialize arrays - only done by root (we pretend this is input
	 data) */
  if(rank==root) {
    for(i=0; i<nrow_a; i++)
      for(j=0;j<ncol_a; j++) 
	a[i][j]=(i+j);
    for(i=0; i<nrow_b; i++)
      for(j=0;j<ncol_b; j++) 
	b[i][j]=(i+j);
  }

  end_init=MPI_Wtime();

  /* Calculate the product in parallel */

  /* Distribute matrix A by rows */

  ierr=MPI_Scatter(*a,nrow_a*ncol_a/size,MPI_DOUBLE,
				   *apart,nrow_a*ncol_a/size,MPI_DOUBLE,
				   root, MPI_COMM_WORLD);

  /* Broadcast matrix B to everybody */

  ierr=MPI_Bcast(*b,nrow_b*ncol_b,MPI_DOUBLE, root, MPI_COMM_WORLD);

  /* Each processor works on its subset of the rows of C. */

  for(i=0; i<nrow_c/size; i++)
    for(j=0;j<ncol_c; j++)
      cpart[i][j]=0.0e0;

  for(i=0; i<nrow_c/size; i++)
    for(k=0; k<ncol_a; k++)
      for(j=0;j<ncol_c; j++) 
	cpart[i][j]+=apart[i][k]*b[k][j];

  /* Collect all the rows back to C on root */

  ierr=MPI_Gather(*cpart,(nrow_c/size)*ncol_c,MPI_DOUBLE,
		  *c,(nrow_c/size)*ncol_c,MPI_DOUBLE,
		  root, MPI_COMM_WORLD);

  /* Root does more I/O and timing, arrays get deallocated. */
  if(rank==root) {

    end_calc=MPI_Wtime();

    printf("c[0][0]=%g exact[0][0]=%g diff=%g\n",
	   c[0][0],
	   exact(0,0,ncol_a),
	   c[0][0]-exact(0,0,ncol_a));

    printf("c[%d][%d]=%g exact[%d][%d]=%g diff=%g\n",
	   nrow_c/size-1,ncol_c-1,c[nrow_c/size-1][ncol_c-1],
	   nrow_c/size-1,ncol_c-1,exact(nrow_c/size-1,ncol_c-1,ncol_a),
	   c[nrow_c/size-1][ncol_c-1]-exact(nrow_c/size-1,ncol_c-1,ncol_a));

    printf("c[%d][%d]=%g exact[%d][%d]=%g diff=%g\n",
	   (nrow_c-1)/2,(ncol_c-1)/2,c[(nrow_c-1)/2][(ncol_c-1)/2],
	   (nrow_c-1)/2,(ncol_c-1)/2,exact((nrow_c-1)/2,(ncol_c-1)/2,ncol_a),
	   c[(nrow_c-1)/2][(ncol_c-1)/2]-exact((nrow_c-1)/2,(ncol_c-1)/2,ncol_a));

    printf("c[%d][%d]=%g exact[%d][%d]=%g diff=%g\n",
	   nrow_c-1,ncol_c-1,c[nrow_c-1][ncol_c-1],
	   nrow_c-1,ncol_c-1,exact(nrow_c-1,ncol_c-1,ncol_a),
	   c[nrow_c-1][ncol_c-1]-exact(nrow_c-1,ncol_c-1,ncol_a));

    deallocate_array(a, nrow_a);
    deallocate_array(c, nrow_c);

  }

  deallocate_array(cpart, nrow_c/size);
  deallocate_array(apart, nrow_a/size);
  deallocate_array(b, nrow_b);

  if(rank==root) {

    end_time=MPI_Wtime();

    printf("Startup time: %18.10lg\n", (double)(end_init-start_time));
    printf("Calculation time: %18.10lg\n", (double)(end_calc-end_init));
    printf("Cleanup time: %18.10lg\n", (double)(end_time-end_calc));
    printf("Total time: %18.10lg\n", (double)(end_time-start_time));

  }

  MPI_Finalize();

}