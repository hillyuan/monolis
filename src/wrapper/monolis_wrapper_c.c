#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "monolis.h"
#include "metis.h"

void monolis_set_method   (MONOLIS* mat, int    flag) {mat->prm.method    = flag;}
void monolis_set_precond  (MONOLIS* mat, int    flag) {mat->prm.precond   = flag;}
void monolis_set_maxiter  (MONOLIS* mat, int    flag) {mat->prm.maxiter   = flag;}
void monolis_set_tolerance(MONOLIS* mat, double flag) {mat->prm.tol       = flag;}
void monolis_show_iterlog (MONOLIS* mat, bool   flag) {mat->prm.show_iterlog = flag;}
void monolis_show_timelog (MONOLIS* mat, bool   flag) {mat->prm.show_timelog = flag;}
void monolis_show_summary (MONOLIS* mat, bool   flag) {mat->prm.show_summary = flag;}

/* initializer */

void monolis_initialize(
  MONOLIS* mat)
{
    // prm
    mat->prm.method = 1;
    mat->prm.precond = 1;
    mat->prm.maxiter = 1000;
    mat->prm.curiter = 0;
    mat->prm.ierr = -1;
    mat->prm.tol = 1.0e-8;
    mat->prm.curresid = 0.0;
    mat->prm.is_scaling = false;
    mat->prm.is_reordering = false;
    mat->prm.is_init_x = true;
    mat->prm.is_debug = false;
    mat->prm.show_iterlog = false;
    mat->prm.show_timelog = false;
    mat->prm.show_summary = true;

    mat->prm.tsol  = 0.0;
    mat->prm.tprep = 0.0;
    mat->prm.tspmv = 0.0;
    mat->prm.tdotp = 0.0;
    mat->prm.tprec = 0.0;
    mat->prm.tcomm = 0.0;

    // mat
    mat->mat.N = 0;
    mat->mat.NP = 0;
    mat->mat.NZ = 0;
    mat->mat.NDOF = 0;

    // comm
    mat->com.myrank = 0;
    mat->com.comm = 0;
    mat->com.commsize = 1;
    mat->com.recv_n_neib = 0;
    mat->com.send_n_neib = 0;
}

void monolis_finalize(
  MONOLIS* mat)
{

}

/* mat clear */

void monolis_clear(
  MONOLIS* mat)
{
  monolis_clear_mat(mat);
  monolis_clear_rhs(mat);
  monolis_clear_solution(mat);
}

void monolis_clear_mat(
  MONOLIS* mat)
{
  int ndof2 = mat->mat.NDOF*mat->mat.NDOF;

  for(int i=0; i<(mat->mat.NZ*ndof2); i++) {
    mat->mat.A[i] = 0.0;
  }
}

void monolis_clear_rhs(
  MONOLIS* mat)
{
  int ndof = mat->mat.NDOF;

  for(int i=0; i<(mat->mat.NP*ndof); i++) {
    mat->mat.B[i] = 0.0;
  }
}

void monolis_clear_solution(
  MONOLIS* mat)
{
  int ndof = mat->mat.NDOF;

  for(int i=0; i<(mat->mat.NP*ndof); i++) {
    mat->mat.X[i] = 0.0;
  }
}

/* mat copy */

void monolis_copy_all(
  MONOLIS* in,
  MONOLIS* out)
{
  int ndof = in->mat.NDOF;
  int np = in->mat.NP;
  int nz = in->mat.NZ;

  monolis_copy_nonzero_pattern(in, out);

  for(int i=0; i<ndof*ndof*nz; i++) {
    out->mat.A[i] = in->mat.A[i];
  }

  for(int i=0; i<ndof*np; i++) {
    out->mat.X[i] = in->mat.X[i];
  }

  for(int i=0; i<ndof*np; i++) {
    out->mat.B[i] = in->mat.B[i];
  }
}

void monolis_copy_nonzero_pattern(
  MONOLIS* in,
  MONOLIS* out)
{
  int n = in->mat.N;
  int np = in->mat.NP;
  int ndof = in->mat.NDOF;
  int nz = in->mat.NZ;

  out->mat.N = n;
  out->mat.NP = np;
  out->mat.NDOF = ndof;
  out->mat.NZ = nz;

  out->mat.A = (double*)calloc(ndof*ndof*nz, sizeof(double));
  for(int i=0; i<ndof*ndof*nz; i++) {
    out->mat.A[i] = 0.0;
  }

  out->mat.X = (double*)calloc(ndof*np, sizeof(double));
  for(int i=0; i<ndof*np; i++) {
    out->mat.X[i] = 0.0;
  }

  out->mat.B = (double*)calloc(ndof*np, sizeof(double));
  for(int i=0; i<ndof*np; i++) {
    out->mat.B[i] = 0.0;
  }

  out->mat.index = (int* )calloc(np+1, sizeof(int));
  for(int i=0; i<np+1; i++) {
    out->mat.index[i] = in->mat.index[i];
  }

  out->mat.item = (int*)calloc(nz, sizeof(int));
  for(int i=0; i<nz; i++) {
    out->mat.item[i] = in->mat.item[i];
  }

  out->mat.indexR = (int*)calloc(np+1, sizeof(int));
  for(int i=0; i<np+1; i++) {
    out->mat.indexR[i] = in->mat.indexR[i];
  }

  out->mat.itemR = (int*)calloc(nz, sizeof(int));
  for(int i=0; i<nz; i++) {
    out->mat.itemR[i] = in->mat.itemR[i];
  }

  out->mat.permR = (int*)calloc(nz, sizeof(int));
  for(int i=0; i<nz; i++) {
    out->mat.permR[i] = in->mat.permR[i];
  }
}

void monolis_get_input_filename(
  )
{
/*
    int commsize;
    int myrank;
    char str_myrank[8];

    commsize = monolis_get_global_commsize();
    myrank = monolis_get_global_myrank();
    snprintf(str_myrank, sizeof(str_myrank), "%d", myrank);

    if( commsize > 1){
      strcpy(INPUT_FILENAME_NODE_PAR, INPUT_FILENAME_NODE);
      strcpy(INPUT_FILENAME_ELEM_PAR, INPUT_FILENAME_ELEM);
    } else {
      strcpy(INPUT_FILENAME_NODE_PAR, INPUT_FILENAME_NODE);
      strcpy(INPUT_FILENAME_ELEM_PAR, INPUT_FILENAME_ELEM);
    }
*/
}

void monolis_convert_mesh_to_connectivity(
  int      nelem,
  int      nbase_func,
  int**    elem,
  idx_t*  conn_index,
  idx_t*  con)
{
  for(int i=0; i<nelem+1; i++) {
    conn_index[i] = i*nbase_func;
  }

  for(int i=0; i<nelem; i++) {
    for(int j=0; j<nbase_func; j++) {
      con[nbase_func*i+j] = elem[i][j];
    }
  }
}

void monolis_convert_connectivity_to_nodal(
  int      nnode,
  int      nelem,
  idx_t*   conn_index,
  idx_t*   con,
  idx_t** index,
  idx_t** item)
{
#ifdef WITH_METIS
  idx_t numflag = 0;
  int ierr = METIS_MeshToNodal(
    &nelem,
    &nnode,
    conn_index,
    con,
    &numflag,
    index,
    item);
#else

#endif
}

void monolis_get_nonzero_pattern_by_nodal(
  MONOLIS* mat,
  int      nnode,
  int      ndof,
  int*     index,
  int*     item)
{
  mat->mat.N = nnode;
  mat->mat.NP = nnode;
  mat->mat.NDOF = ndof;
  mat->mat.X = (double*)calloc(ndof*nnode, sizeof(double));
  mat->mat.B = (double*)calloc(ndof*nnode, sizeof(double));

  mat->mat.index = (int* )calloc(nnode+1, sizeof(int));
  for(int i=1; i<nnode+1; i++) {
    mat->mat.index[i] = index[i] + i;
  }

  int nz = mat->mat.index[nnode];
  mat->mat.NZ = nz;
  mat->mat.A = (double*)calloc(ndof*ndof*nz, sizeof(double));
  mat->mat.item = (int*)calloc(nz, sizeof(int));

  for(int i=0; i<nnode; i++) {
    int jS = mat->mat.index[i];
    int jE = mat->mat.index[i+1];
    mat->mat.item[jS] = i+1;
    for(int j=jS+1; j<jE+1; j++){
      mat->mat.item[j] = item[j-i-1] + 1;
    }
    monolis_qsort_int(
      &(mat->mat.item[jS]),
      1,
      jE-jS);
  }

  mat->mat.indexR = (int*)calloc(nnode+1, sizeof(int));
  mat->mat.itemR = (int*)calloc(nz, sizeof(int));
  mat->mat.permR = (int*)calloc(nz, sizeof(int));
  monolis_get_CRR_format(
    nnode,
    nz,
    mat->mat.index,
    mat->mat.item,
    mat->mat.indexR,
    mat->mat.itemR,
    mat->mat.permR);
}

void monolis_get_nonzero_pattern(
  MONOLIS* mat,
  int      nnode,
  int      nbase_func,
  int      ndof,
  int      nelem,
  int**    elem)
{
  idx_t* conn_index;
  idx_t* con;
  idx_t* index;
  idx_t* item;

  conn_index = (idx_t*)calloc(nelem+1, sizeof(idx_t));
  con = (idx_t*)calloc(nelem*nbase_func, sizeof(idx_t));

  monolis_convert_mesh_to_connectivity(
    nelem,
    nbase_func,
    elem,
    conn_index,
    con);

  monolis_convert_connectivity_to_nodal(
    nnode,
    nelem,
    conn_index,
    con,
    &index,
    &item);

  monolis_get_nonzero_pattern_by_nodal(
    mat,
    nnode,
    ndof,
    index,
    item);

  free(conn_index);
  free(con);
  free(index);
  free(item);
}

void monolis_add_scalar_to_sparse_matrix(
  MONOLIS* mat,
  double   val,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j)
{
  int nnode = mat->mat.N;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.NZ;

  monolis_add_scalar_to_sparse_matrix_c_main(
    nnode,
    nz,
    ndof,
    mat->mat.index,
    mat->mat.item,
    mat->mat.A,
    i,
    j,
    submat_i,
    submat_j,
    val);
}

void monolis_set_Dirichlet_bc(
  MONOLIS* mat,
  double*  b,
  int      node_id,
  int      ndof_bc,
  double   val)
{
  int nnode = mat->mat.N;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.NZ;

  monolis_set_Dirichlet_bc_c_main(
    nnode,
    nz,
    ndof,
    mat->mat.index,
    mat->mat.item,
    mat->mat.indexR,
    mat->mat.itemR,
    mat->mat.permR,
    mat->mat.A,
    b,
    node_id,
    ndof_bc,
    val);
}

void monolis_solve(
  MONOLIS* mat,
  double*  b,
  double*  x)
{
  int n = mat->mat.N;
  int np = mat->mat.N;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.NZ;
  int iterlog = 0;
  int timelog = 0;
  int summary = 0;

  if(mat->prm.show_iterlog) iterlog = 1;
  if(mat->prm.show_timelog) timelog = 1;
  if(mat->prm.show_summary) summary = 1;

  monolis_solve_c_main(
    n,
    np,
    nz,
    ndof,
    mat->mat.A,
    x,
    b,
    mat->mat.index,
    mat->mat.item,
    mat->com.myrank,
    mat->com.comm,
    mat->com.commsize,
    mat->com.recv_n_neib,
    mat->com.send_n_neib,
    mat->prm.method,
    mat->prm.precond,
    mat->prm.maxiter,
    mat->prm.tol,
    iterlog,
    timelog,
    summary);
}
