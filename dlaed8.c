#include "rb_lapack.h"

extern VOID dlaed8_(integer *icompq, integer *k, integer *n, integer *qsiz, doublereal *d, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, integer *cutpnt, doublereal *z, doublereal *dlamda, doublereal *q2, integer *ldq2, doublereal *w, integer *perm, integer *givptr, integer *givcol, doublereal *givnum, integer *indxp, integer *indx, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlaed8(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_icompq;
  integer icompq; 
  VALUE rblapack_qsiz;
  integer qsiz; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_q;
  doublereal *q; 
  VALUE rblapack_ldq;
  integer ldq; 
  VALUE rblapack_indxq;
  integer *indxq; 
  VALUE rblapack_rho;
  doublereal rho; 
  VALUE rblapack_cutpnt;
  integer cutpnt; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_ldq2;
  integer ldq2; 
  VALUE rblapack_k;
  integer k; 
  VALUE rblapack_dlamda;
  doublereal *dlamda; 
  VALUE rblapack_q2;
  doublereal *q2; 
  VALUE rblapack_w;
  doublereal *w; 
  VALUE rblapack_perm;
  integer *perm; 
  VALUE rblapack_givptr;
  integer givptr; 
  VALUE rblapack_givcol;
  integer *givcol; 
  VALUE rblapack_givnum;
  doublereal *givnum; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  VALUE rblapack_q_out__;
  doublereal *q_out__;
  integer *indxp;
  integer *indx;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  k, dlamda, q2, w, perm, givptr, givcol, givnum, info, d, q, rho = NumRu::Lapack.dlaed8( icompq, qsiz, d, q, ldq, indxq, rho, cutpnt, z, ldq2, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO, CUTPNT, Z, DLAMDA, Q2, LDQ2, W, PERM, GIVPTR, GIVCOL, GIVNUM, INDXP, INDX, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLAED8 merges the two sets of eigenvalues together into a single\n*  sorted set.  Then it tries to deflate the size of the problem.\n*  There are two ways in which deflation can occur:  when two or more\n*  eigenvalues are close together or if there is a tiny element in the\n*  Z vector.  For each such occurrence the order of the related secular\n*  equation problem is reduced by one.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ  (input) INTEGER\n*          = 0:  Compute eigenvalues only.\n*          = 1:  Compute eigenvectors of original dense symmetric matrix\n*                also.  On entry, Q contains the orthogonal matrix used\n*                to reduce the original matrix to tridiagonal form.\n*\n*  K      (output) INTEGER\n*         The number of non-deflated eigenvalues, and the order of the\n*         related secular equation.\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  QSIZ   (input) INTEGER\n*         The dimension of the orthogonal matrix used to reduce\n*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension (N)\n*         On entry, the eigenvalues of the two submatrices to be\n*         combined.  On exit, the trailing (N-K) updated eigenvalues\n*         (those which were deflated) sorted into increasing order.\n*\n*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ,N)\n*         If ICOMPQ = 0, Q is not referenced.  Otherwise,\n*         on entry, Q contains the eigenvectors of the partially solved\n*         system which has been previously updated in matrix\n*         multiplies with other partially solved eigensystems.\n*         On exit, Q contains the trailing (N-K) updated eigenvectors\n*         (those which were deflated) in its last N-K columns.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  INDXQ  (input) INTEGER array, dimension (N)\n*         The permutation which separately sorts the two sub-problems\n*         in D into ascending order.  Note that elements in the second\n*         half of this permutation must first have CUTPNT added to\n*         their values in order to be accurate.\n*\n*  RHO    (input/output) DOUBLE PRECISION\n*         On entry, the off-diagonal element associated with the rank-1\n*         cut which originally split the two submatrices which are now\n*         being recombined.\n*         On exit, RHO has been modified to the value required by\n*         DLAED3.\n*\n*  CUTPNT (input) INTEGER\n*         The location of the last eigenvalue in the leading\n*         sub-matrix.  min(1,N) <= CUTPNT <= N.\n*\n*  Z      (input) DOUBLE PRECISION array, dimension (N)\n*         On entry, Z contains the updating vector (the last row of\n*         the first sub-eigenvector matrix and the first row of the\n*         second sub-eigenvector matrix).\n*         On exit, the contents of Z are destroyed by the updating\n*         process.\n*\n*  DLAMDA (output) DOUBLE PRECISION array, dimension (N)\n*         A copy of the first K eigenvalues which will be used by\n*         DLAED3 to form the secular equation.\n*\n*  Q2     (output) DOUBLE PRECISION array, dimension (LDQ2,N)\n*         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,\n*         a copy of the first K eigenvectors which will be used by\n*         DLAED7 in a matrix multiply (DGEMM) to update the new\n*         eigenvectors.\n*\n*  LDQ2   (input) INTEGER\n*         The leading dimension of the array Q2.  LDQ2 >= max(1,N).\n*\n*  W      (output) DOUBLE PRECISION array, dimension (N)\n*         The first k values of the final deflation-altered z-vector and\n*         will be passed to DLAED3.\n*\n*  PERM   (output) INTEGER array, dimension (N)\n*         The permutations (from deflation and sorting) to be applied\n*         to each eigenblock.\n*\n*  GIVPTR (output) INTEGER\n*         The number of Givens rotations which took place in this\n*         subproblem.\n*\n*  GIVCOL (output) INTEGER array, dimension (2, N)\n*         Each pair of numbers indicates a pair of columns to take place\n*         in a Givens rotation.\n*\n*  GIVNUM (output) DOUBLE PRECISION array, dimension (2, N)\n*         Each number indicates the S value to be used in the\n*         corresponding Givens rotation.\n*\n*  INDXP  (workspace) INTEGER array, dimension (N)\n*         The permutation used to place deflated values of D at the end\n*         of the array.  INDXP(1:K) points to the nondeflated D-values\n*         and INDXP(K+1:N) points to the deflated eigenvalues.\n*\n*  INDX   (workspace) INTEGER array, dimension (N)\n*         The permutation used to sort the contents of D into ascending\n*         order.\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  k, dlamda, q2, w, perm, givptr, givcol, givnum, info, d, q, rho = NumRu::Lapack.dlaed8( icompq, qsiz, d, q, ldq, indxq, rho, cutpnt, z, ldq2, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_icompq = argv[0];
  rblapack_qsiz = argv[1];
  rblapack_d = argv[2];
  rblapack_q = argv[3];
  rblapack_ldq = argv[4];
  rblapack_indxq = argv[5];
  rblapack_rho = argv[6];
  rblapack_cutpnt = argv[7];
  rblapack_z = argv[8];
  rblapack_ldq2 = argv[9];
  if (rb_options != Qnil) {
  }

  qsiz = NUM2INT(rblapack_qsiz);
  cutpnt = NUM2INT(rblapack_cutpnt);
  if (!NA_IsNArray(rblapack_indxq))
    rb_raise(rb_eArgError, "indxq (6th argument) must be NArray");
  if (NA_RANK(rblapack_indxq) != 1)
    rb_raise(rb_eArgError, "rank of indxq (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_indxq);
  if (NA_TYPE(rblapack_indxq) != NA_LINT)
    rblapack_indxq = na_change_type(rblapack_indxq, NA_LINT);
  indxq = NA_PTR_TYPE(rblapack_indxq, integer*);
  rho = NUM2DBL(rblapack_rho);
  ldq = NUM2INT(rblapack_ldq);
  icompq = NUM2INT(rblapack_icompq);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of indxq");
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of indxq");
  if (NA_TYPE(rblapack_z) != NA_DFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (4th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != (icompq==0 ? 0 : n))
    rb_raise(rb_eRuntimeError, "shape 1 of q must be %d", icompq==0 ? 0 : n);
  if (NA_SHAPE0(rblapack_q) != (icompq==0 ? 0 : ldq))
    rb_raise(rb_eRuntimeError, "shape 0 of q must be %d", icompq==0 ? 0 : ldq);
  if (NA_TYPE(rblapack_q) != NA_DFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rblapack_q, doublereal*);
  ldq2 = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rblapack_dlamda = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dlamda = NA_PTR_TYPE(rblapack_dlamda, doublereal*);
  {
    int shape[2];
    shape[0] = icompq==0 ? 0 : ldq2;
    shape[1] = icompq==0 ? 0 : n;
    rblapack_q2 = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q2 = NA_PTR_TYPE(rblapack_q2, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_perm = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  perm = NA_PTR_TYPE(rblapack_perm, integer*);
  {
    int shape[2];
    shape[0] = 2;
    shape[1] = n;
    rblapack_givcol = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  givcol = NA_PTR_TYPE(rblapack_givcol, integer*);
  {
    int shape[2];
    shape[0] = 2;
    shape[1] = n;
    rblapack_givnum = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rblapack_givnum, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = icompq==0 ? 0 : ldq;
    shape[1] = icompq==0 ? 0 : n;
    rblapack_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  indxp = ALLOC_N(integer, (n));
  indx = ALLOC_N(integer, (n));

  dlaed8_(&icompq, &k, &n, &qsiz, d, q, &ldq, indxq, &rho, &cutpnt, z, dlamda, q2, &ldq2, w, perm, &givptr, givcol, givnum, indxp, indx, &info);

  free(indxp);
  free(indx);
  rblapack_k = INT2NUM(k);
  rblapack_givptr = INT2NUM(givptr);
  rblapack_info = INT2NUM(info);
  rblapack_rho = rb_float_new((double)rho);
  return rb_ary_new3(12, rblapack_k, rblapack_dlamda, rblapack_q2, rblapack_w, rblapack_perm, rblapack_givptr, rblapack_givcol, rblapack_givnum, rblapack_info, rblapack_d, rblapack_q, rblapack_rho);
}

void
init_lapack_dlaed8(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed8", rblapack_dlaed8, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
