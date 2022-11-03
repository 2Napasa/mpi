/* Algorithm:
 *    1.  Each process calculates "its" interval of
 *        integration.
 *    2.  Each process outs the integral of f(x)
 *        over its interval using the trapezoidalal rule.
 *    3a. Each process != 0 sends its integral to 0.
 *    3b. Process 0 sums the calculations received from
 *        the individual processes and prints the result.
*/
#include <stdio.h>
#include <mpi.h>
#include <math.h>
double trapezoidal(double left, double right, int count, double dx);    
double f(double x); 

int main(void) {
   int my_rank, comm_sz, n=5000, local_n;   
   double a=1, b=10, h, local_a, local_b;
   double local_int, total_int;
   int source; 
   
   MPI_Init(NULL, NULL);                        /* Let the system do what it needs to start up MPI */
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);     /* Get my process rank */
   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);     /* Find out how many processes are being used */
   
   h = (b-a)/n;                                 /* h is the same for all processes */
   local_n = n/comm_sz;                         /* So is the number of trapezoidals  */

   /* Length of each process' interval of
    * integration = local_n*h.  So my interval
    * starts at: */
   local_a = a + my_rank*local_n*h;
   local_b = local_a + local_n*h;
   local_int = trapezoidal(local_a, local_b, local_n, h);

   /* Add up the integrals calculated by each process */
   if (my_rank != 0)
      MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
   else {
      total_int = local_int;
      for (source = 1; source < comm_sz; source++) {
         MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         total_int += local_int;
      }
   }

   if (my_rank == 0) {
      printf(" = %.6f\n", total_int);
   }
   MPI_Finalize();
   return 0;
} 

double trapezoidal(double left, double right, int count , double dx) {
   double out, x; 
   int i;
   out = (f(left) + f(right))/2.0;
   for (i = 1; i <= count-1; i++) {
      x = left + i*dx;
      out += f(x);
   }
   out = out*dx;
   return out;
} 

double f(double x) {
   return x*x*sin(x)*sin(x);
} 