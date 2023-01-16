#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "monolis.h"

/*
We have investigated the problem further, and we have
been able to reproduce it and obtain an erroneous
solution with an even smaller, 2x2, matrix:
    [2 1]
    [1 3]
and a right-hand side vector with all ones (1,1)
The correct solution is the vector (0.4,0.2), in both solves.

mpirun -n 1 ./main
*/

int main(int argc, char *args[]) {
  MONOLIS monolis;
  int nz;
  double res;
  double A[4], b[2], sol[2];
  int index[3];
  int jcol[4];
  int gindex[2];
  
  monolis_global_initialize();

  monolis_prm_initialize(&monolis);
  monolis_com_initialize(&monolis);
  monolis_mat_initialize(&monolis);
  
  assert(monolis.com.commsize==1);

  index[0] = 0; index[1]=2; index[2]=4;  nz=4;
  gindex[0]=1;  gindex[1]=2;
  jcol[0]=0; jcol[1]=1;  jcol[2]=0; jcol[3]=1;
  A[0]=2.0; A[1]=1.0;  A[2]=1.0; A[3]=3.0;

  monolis_set_matrix_BCSR(&monolis,2,2,1,nz,A,index,jcol);
  monolis_com_get_comm_table(&monolis,2,2,gindex);

  sol[0]= 0.0;  sol[1]=0.0;
  b[0]=1.0;  b[1]=1.0;
  monolis_param_set_show_iterlog(&monolis,true);
  monolis_solve(&monolis,b,sol);
  
  monolis_get_converge_residual(&monolis, &res);
  printf("* monolis_get_converge_residual     %e\n", res);

  printf("%f\n", sol[0]);
  printf("%f\n", sol[1]);

  monolis_finalize(&monolis);
  monolis_global_finalize();
}
