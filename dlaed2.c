#include "rb_lapack.h"

extern VOID dlaed2_(integer *k, integer *n, integer *n1, doublereal *d, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, doublereal *z, doublereal *dlamda, doublereal *w, doublereal *q2, integer *indx, integer *indxc, integer *indxp, integer *coltyp, integer *info);

static VALUE
rb_dlaed2(int argc, VALUE *argv, VALUE self){
  VALUE rb_n1;
  integer n1; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_rho;
  doublereal rho; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_k;
  integer k; 
  VALUE rb_dlamda;
  doublereal *dlamda; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_q2;
  doublereal *q2; 
  VALUE rb_indxc;
  integer *indxc; 
  VALUE rb_coltyp;
  integer *coltyp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_q_out__;
  doublereal *q_out__;
  VALUE rb_indxq_out__;
  integer *indxq_out__;
  integer *indx;
  integer *indxp;

  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, dlamda, w, q2, indxc, coltyp, info, d, q, indxq, rho = NumRu::Lapack.dlaed2( n1, d, q, indxq, rho, z)\n    or\n  NumRu::Lapack.dlaed2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMDA, W, Q2, INDX, INDXC, INDXP, COLTYP, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLAED2 merges the two sets of eigenvalues together into a single\n*  sorted set.  Then it tries to deflate the size of the problem.\n*  There are two ways in which deflation can occur:  when two or more\n*  eigenvalues are close together or if there is a tiny entry in the\n*  Z vector.  For each such occurrence the order of the related secular\n*  equation problem is reduced by one.\n*\n\n*  Arguments\n*  =========\n*\n*  K      (output) INTEGER\n*         The number of non-deflated eigenvalues, and the order of the\n*         related secular equation. 0 <= K <=N.\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  N1     (input) INTEGER\n*         The location of the last eigenvalue in the leading sub-matrix.\n*         min(1,N) <= N1 <= N/2.\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension (N)\n*         On entry, D contains the eigenvalues of the two submatrices to\n*         be combined.\n*         On exit, D contains the trailing (N-K) updated eigenvalues\n*         (those which were deflated) sorted into increasing order.\n*\n*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)\n*         On entry, Q contains the eigenvectors of two submatrices in\n*         the two square blocks with corners at (1,1), (N1,N1)\n*         and (N1+1, N1+1), (N,N).\n*         On exit, Q contains the trailing (N-K) updated eigenvectors\n*         (those which were deflated) in its last N-K columns.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  INDXQ  (input/output) INTEGER array, dimension (N)\n*         The permutation which separately sorts the two sub-problems\n*         in D into ascending order.  Note that elements in the second\n*         half of this permutation must first have N1 added to their\n*         values. Destroyed on exit.\n*\n*  RHO    (input/output) DOUBLE PRECISION\n*         On entry, the off-diagonal element associated with the rank-1\n*         cut which originally split the two submatrices which are now\n*         being recombined.\n*         On exit, RHO has been modified to the value required by\n*         DLAED3.\n*\n*  Z      (input) DOUBLE PRECISION array, dimension (N)\n*         On entry, Z contains the updating vector (the last\n*         row of the first sub-eigenvector matrix and the first row of\n*         the second sub-eigenvector matrix).\n*         On exit, the contents of Z have been destroyed by the updating\n*         process.\n*\n*  DLAMDA (output) DOUBLE PRECISION array, dimension (N)\n*         A copy of the first K eigenvalues which will be used by\n*         DLAED3 to form the secular equation.\n*\n*  W      (output) DOUBLE PRECISION array, dimension (N)\n*         The first k values of the final deflation-altered z-vector\n*         which will be passed to DLAED3.\n*\n*  Q2     (output) DOUBLE PRECISION array, dimension (N1**2+(N-N1)**2)\n*         A copy of the first K eigenvectors which will be used by\n*         DLAED3 in a matrix multiply (DGEMM) to solve for the new\n*         eigenvectors.\n*\n*  INDX   (workspace) INTEGER array, dimension (N)\n*         The permutation used to sort the contents of DLAMDA into\n*         ascending order.\n*\n*  INDXC  (output) INTEGER array, dimension (N)\n*         The permutation used to arrange the columns of the deflated\n*         Q matrix into three groups:  the first group contains non-zero\n*         elements only at and above N1, the second contains\n*         non-zero elements only below N1, and the third is dense.\n*\n*  INDXP  (workspace) INTEGER array, dimension (N)\n*         The permutation used to place deflated values of D at the end\n*         of the array.  INDXP(1:K) points to the nondeflated D-values\n*         and INDXP(K+1:N) points to the deflated eigenvalues.\n*\n*  COLTYP (workspace/output) INTEGER array, dimension (N)\n*         During execution, a label which will indicate which of the\n*         following types a column in the Q2 matrix is:\n*         1 : non-zero in the upper half only;\n*         2 : dense;\n*         3 : non-zero in the lower half only;\n*         4 : deflated.\n*         On exit, COLTYP(i) is the number of columns of type i,\n*         for i=1 to 4 only.\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*  Modified by Francoise Tisseur, University of Tennessee.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_n1 = argv[0];
  rb_d = argv[1];
  rb_q = argv[2];
  rb_indxq = argv[3];
  rb_rho = argv[4];
  rb_z = argv[5];

  if (!NA_IsNArray(rb_indxq))
    rb_raise(rb_eArgError, "indxq (4th argument) must be NArray");
  if (NA_RANK(rb_indxq) != 1)
    rb_raise(rb_eArgError, "rank of indxq (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_indxq);
  if (NA_TYPE(rb_indxq) != NA_LINT)
    rb_indxq = na_change_type(rb_indxq, NA_LINT);
  indxq = NA_PTR_TYPE(rb_indxq, integer*);
  rho = NUM2DBL(rb_rho);
  n1 = NUM2INT(rb_n1);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (6th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (3th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of indxq");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_dlamda = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dlamda = NA_PTR_TYPE(rb_dlamda, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[1];
    shape[0] = pow(n1,2)+pow(n-n1,2);
    rb_q2 = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  q2 = NA_PTR_TYPE(rb_q2, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_indxc = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indxc = NA_PTR_TYPE(rb_indxc, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_coltyp = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  coltyp = NA_PTR_TYPE(rb_coltyp, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_indxq_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indxq_out__ = NA_PTR_TYPE(rb_indxq_out__, integer*);
  MEMCPY(indxq_out__, indxq, integer, NA_TOTAL(rb_indxq));
  rb_indxq = rb_indxq_out__;
  indxq = indxq_out__;
  indx = ALLOC_N(integer, (n));
  indxp = ALLOC_N(integer, (n));

  dlaed2_(&k, &n, &n1, d, q, &ldq, indxq, &rho, z, dlamda, w, q2, indx, indxc, indxp, coltyp, &info);

  free(indx);
  free(indxp);
  rb_k = INT2NUM(k);
  rb_info = INT2NUM(info);
  rb_rho = rb_float_new((double)rho);
  return rb_ary_new3(11, rb_k, rb_dlamda, rb_w, rb_q2, rb_indxc, rb_coltyp, rb_info, rb_d, rb_q, rb_indxq, rb_rho);
}

void
init_lapack_dlaed2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed2", rb_dlaed2, -1);
}
