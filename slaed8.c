#include "rb_lapack.h"

static VALUE
rb_slaed8(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_d;
  real *d; 
  VALUE rb_q;
  real *q; 
  VALUE rb_ldq;
  integer ldq; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_cutpnt;
  integer cutpnt; 
  VALUE rb_z;
  real *z; 
  VALUE rb_ldq2;
  integer ldq2; 
  VALUE rb_k;
  integer k; 
  VALUE rb_dlamda;
  real *dlamda; 
  VALUE rb_q2;
  real *q2; 
  VALUE rb_w;
  real *w; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  real *givnum; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_q_out__;
  real *q_out__;
  integer *indxp;
  integer *indx;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, dlamda, q2, w, perm, givptr, givcol, givnum, info, d, q, rho = NumRu::Lapack.slaed8( icompq, qsiz, d, q, ldq, indxq, rho, cutpnt, z, ldq2)\n    or\n  NumRu::Lapack.slaed8  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO, CUTPNT, Z, DLAMDA, Q2, LDQ2, W, PERM, GIVPTR, GIVCOL, GIVNUM, INDXP, INDX, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAED8 merges the two sets of eigenvalues together into a single\n*  sorted set.  Then it tries to deflate the size of the problem.\n*  There are two ways in which deflation can occur:  when two or more\n*  eigenvalues are close together or if there is a tiny element in the\n*  Z vector.  For each such occurrence the order of the related secular\n*  equation problem is reduced by one.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ  (input) INTEGER\n*          = 0:  Compute eigenvalues only.\n*          = 1:  Compute eigenvectors of original dense symmetric matrix\n*                also.  On entry, Q contains the orthogonal matrix used\n*                to reduce the original matrix to tridiagonal form.\n*\n*  K      (output) INTEGER\n*         The number of non-deflated eigenvalues, and the order of the\n*         related secular equation.\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  QSIZ   (input) INTEGER\n*         The dimension of the orthogonal matrix used to reduce\n*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.\n*\n*  D      (input/output) REAL array, dimension (N)\n*         On entry, the eigenvalues of the two submatrices to be\n*         combined.  On exit, the trailing (N-K) updated eigenvalues\n*         (those which were deflated) sorted into increasing order.\n*\n*  Q      (input/output) REAL array, dimension (LDQ,N)\n*         If ICOMPQ = 0, Q is not referenced.  Otherwise,\n*         on entry, Q contains the eigenvectors of the partially solved\n*         system which has been previously updated in matrix\n*         multiplies with other partially solved eigensystems.\n*         On exit, Q contains the trailing (N-K) updated eigenvectors\n*         (those which were deflated) in its last N-K columns.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  INDXQ  (input) INTEGER array, dimension (N)\n*         The permutation which separately sorts the two sub-problems\n*         in D into ascending order.  Note that elements in the second\n*         half of this permutation must first have CUTPNT added to\n*         their values in order to be accurate.\n*\n*  RHO    (input/output) REAL\n*         On entry, the off-diagonal element associated with the rank-1\n*         cut which originally split the two submatrices which are now\n*         being recombined.\n*         On exit, RHO has been modified to the value required by\n*         SLAED3.\n*\n*  CUTPNT (input) INTEGER\n*         The location of the last eigenvalue in the leading\n*         sub-matrix.  min(1,N) <= CUTPNT <= N.\n*\n*  Z      (input) REAL array, dimension (N)\n*         On entry, Z contains the updating vector (the last row of\n*         the first sub-eigenvector matrix and the first row of the\n*         second sub-eigenvector matrix).\n*         On exit, the contents of Z are destroyed by the updating\n*         process.\n*\n*  DLAMDA (output) REAL array, dimension (N)\n*         A copy of the first K eigenvalues which will be used by\n*         SLAED3 to form the secular equation.\n*\n*  Q2     (output) REAL array, dimension (LDQ2,N)\n*         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,\n*         a copy of the first K eigenvectors which will be used by\n*         SLAED7 in a matrix multiply (SGEMM) to update the new\n*         eigenvectors.\n*\n*  LDQ2   (input) INTEGER\n*         The leading dimension of the array Q2.  LDQ2 >= max(1,N).\n*\n*  W      (output) REAL array, dimension (N)\n*         The first k values of the final deflation-altered z-vector and\n*         will be passed to SLAED3.\n*\n*  PERM   (output) INTEGER array, dimension (N)\n*         The permutations (from deflation and sorting) to be applied\n*         to each eigenblock.\n*\n*  GIVPTR (output) INTEGER\n*         The number of Givens rotations which took place in this\n*         subproblem.\n*\n*  GIVCOL (output) INTEGER array, dimension (2, N)\n*         Each pair of numbers indicates a pair of columns to take place\n*         in a Givens rotation.\n*\n*  GIVNUM (output) REAL array, dimension (2, N)\n*         Each number indicates the S value to be used in the\n*         corresponding Givens rotation.\n*\n*  INDXP  (workspace) INTEGER array, dimension (N)\n*         The permutation used to place deflated values of D at the end\n*         of the array.  INDXP(1:K) points to the nondeflated D-values\n*         and INDXP(K+1:N) points to the deflated eigenvalues.\n*\n*  INDX   (workspace) INTEGER array, dimension (N)\n*         The permutation used to sort the contents of D into ascending\n*         order.\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_icompq = argv[0];
  rb_qsiz = argv[1];
  rb_d = argv[2];
  rb_q = argv[3];
  rb_ldq = argv[4];
  rb_indxq = argv[5];
  rb_rho = argv[6];
  rb_cutpnt = argv[7];
  rb_z = argv[8];
  rb_ldq2 = argv[9];

  icompq = NUM2INT(rb_icompq);
  qsiz = NUM2INT(rb_qsiz);
  ldq = NUM2INT(rb_ldq);
  rho = (real)NUM2DBL(rb_rho);
  cutpnt = NUM2INT(rb_cutpnt);
  ldq2 = NUM2INT(rb_ldq2);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (4th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (4th argument) must be %d", 2);
  if (NA_SHAPE0(rb_q) != (icompq==0 ? 0 : ldq))
    rb_raise(rb_eRuntimeError, "shape 0 of q must be %d", icompq==0 ? 0 : ldq);
  if (NA_SHAPE1(rb_q) != (icompq==0 ? 0 : n))
    rb_raise(rb_eRuntimeError, "shape 1 of q must be %d", icompq==0 ? 0 : n);
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
  if (!NA_IsNArray(rb_indxq))
    rb_raise(rb_eArgError, "indxq (6th argument) must be NArray");
  if (NA_RANK(rb_indxq) != 1)
    rb_raise(rb_eArgError, "rank of indxq (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_indxq) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of indxq must be the same as shape 0 of d");
  if (NA_TYPE(rb_indxq) != NA_LINT)
    rb_indxq = na_change_type(rb_indxq, NA_LINT);
  indxq = NA_PTR_TYPE(rb_indxq, integer*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of d");
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_dlamda = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dlamda = NA_PTR_TYPE(rb_dlamda, real*);
  ldq2 = MAX(1,n);
  {
    int shape[2];
    shape[0] = icompq==0 ? 0 : ldq2;
    shape[1] = icompq==0 ? 0 : n;
    rb_q2 = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q2 = NA_PTR_TYPE(rb_q2, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_perm = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  perm = NA_PTR_TYPE(rb_perm, integer*);
  {
    int shape[2];
    shape[0] = DIM_LEN(2);
    shape[1] = DIM_LEN(n);
    rb_givcol = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  {
    int shape[2];
    shape[0] = 2;
    shape[1] = n;
    rb_givnum = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rb_givnum, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = icompq==0 ? 0 : ldq;
    shape[1] = icompq==0 ? 0 : n;
    rb_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  indxp = ALLOC_N(integer, (n));
  indx = ALLOC_N(integer, (n));

  slaed8_(&icompq, &k, &n, &qsiz, d, q, &ldq, indxq, &rho, &cutpnt, z, dlamda, q2, &ldq2, w, perm, &givptr, givcol, givnum, indxp, indx, &info);

  free(indxp);
  free(indx);
  rb_k = INT2NUM(k);
  rb_givptr = INT2NUM(givptr);
  rb_info = INT2NUM(info);
  rb_rho = rb_float_new((double)rho);
  return rb_ary_new3(12, rb_k, rb_dlamda, rb_q2, rb_w, rb_perm, rb_givptr, rb_givcol, rb_givnum, rb_info, rb_d, rb_q, rb_rho);
}

void
init_lapack_slaed8(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed8", rb_slaed8, -1);
}
