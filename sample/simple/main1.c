#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "monolis.h"

/*
We have investigated the problem further, and we have
been able to reproduce it and obtain an erroneous
solution with an even smaller, 2x2, matrix:
    [2 1]   <---- in cpu0
    [1 3]   <---- in cpu1
and a right-hand side vector with all ones (1,1)
The correct solution is the vector (0.4,0.2), in both solves.

mpirun -n 2 ./main
*/

int main(int argc, char *args[]) {
  MONOLIS monolis;
  int nz;
  double res;
  double A[2], b[2], sol[2];
  int index[2];
  int jcol[2];
  int gindex[2];
  
  monolis_global_initialize();

  monolis_prm_initialize(&monolis);
  monolis_com_initialize(&monolis);
  monolis_mat_initialize(&monolis);
  
  assert(monolis.com.commsize==2);

  index[0] = 0; index[1]=2; nz=2;

  if( monolis.com.myrank==0 ) {
    gindex[0]=1;  gindex[1]=2;
    jcol[0]=0; jcol[1]=1;
    A[0]=2.0; A[1]=1.0;
  } else {
    gindex[0]=2;  gindex[1]=1;
    jcol[0]=0; jcol[1]=1;
    A[0]=3.0; A[1]=1.0;
  }

  monolis_set_matrix_BCSR(&monolis,1,2,1,nz,A,index,jcol);
  monolis_com_get_comm_table(&monolis,1,2,gindex);

  b[0]=1.0; b[1]=1.0;
  monolis_solve(&monolis,b,sol);
  
  monolis_get_converge_residual(&monolis, &res);
  printf("* monolis_get_converge_residual     %e\n", res);

  printf("%f\n", sol[0]);
  printf("%f\n", sol[1]);

  monolis_finalize(&monolis);
  monolis_global_finalize();
}
