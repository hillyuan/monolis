#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "monolis.h"

/*
We have investigated the problem further, and we have
been able to reproduce it and obtain an erroneous
solution with an even smaller, 2x2, matrix:
    [1 2]
    [2 3]
and a right-hand side vector with all ones (1,1)
The correct solution is the vector (-1,1), in both solves.

mpirun -n 2 ./main
*/

int main(int argc, char *args[]) {
  MONOLIS monolis;
  const char* dir_name;
  int iter,nz;
  double time[7], res;
  double A[4], b[2], sol[2];
  int index[3];
  int jcol[4];
  int gindex[2];
  
  monolis_global_initialize();

  monolis_prm_initialize(&monolis);
  monolis_com_initialize(&monolis);
  monolis_mat_initialize(&monolis);
  
  assert(monolis.com.commsize==2);

  index[0] = 0; index[1]=2; index[2]=4; nz=4;
  if( monolis.com.commsize==0 ) {
    gindex[0]=1;  gindex[2]=2;
	jcol[0]=0; jcol[1]=1;  jcol[2]=0; jcol[3]=1;
	A[0]=1.0; A[1]=2.0;  A[2]=3.0;  A[3]=4.0;
  } else {
	gindex[1]=2;  gindex[2]=1; 
    jcol[0]=1; jcol[1]=0;  jcol[2]=1; jcol[3]=0;
    A[0]=2.0; A[1]=1.0;  A[2]=4.0;  A[3]=3.0;
  }
  monolis_set_matrix_BCSR(&monolis,1,2,1,nz,A,index,jcol);
  monolis_com_get_comm_table(&monolis,1,2,gindex);

  b[0]=1.0; b[1]=1.0;
  monolis_solve(&monolis,b,sol);

  monolis_get_time_solver            (&monolis, &time[0]);
  monolis_get_time_preparing         (&monolis, &time[1]);
  monolis_get_time_spmv              (&monolis, &time[2]);
  monolis_get_time_inner_product     (&monolis, &time[3]);
  monolis_get_time_precondition      (&monolis, &time[4]);
  monolis_get_time_comm_inner_product(&monolis, &time[5]);
  monolis_get_time_comm_spmv         (&monolis, &time[6]);
  monolis_get_converge_iter          (&monolis, &iter);
  monolis_get_converge_residual      (&monolis, &res);


  printf("* monolis_get_time_solver           %e\n", time[0]);
  printf("* monolis_get_time_preparing        %e\n", time[1]);
  printf("* monolis_get_time_spmv             %e\n", time[2]);
  printf("* monolis_get_time_dot_product      %e\n", time[3]);
  printf("* monolis_get_time_precondition     %e\n", time[4]);
  printf("* monolis_get_time_comm_dot_product %e\n", time[5]);
  printf("* monolis_get_time_comm_spmv        %e\n", time[6]);
  printf("* monolis_get_converge_iter         %d\n", iter);
  printf("* monolis_get_converge_residual     %e\n", res);

  printf("* monolis_finalize\n");
  monolis_finalize(&monolis);

  printf("* monolis_global_finalize\n");
  monolis_global_finalize();
}
