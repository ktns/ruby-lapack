#include "rb_lapack.h"

static VALUE
rb_zlaed8(int argc, VALUE *argv, VALUE self){
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_q;
  doublecomplex *q; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_rho;
  doublereal rho; 
  VALUE rb_cutpnt;
  integer cutpnt; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_k;
  integer k; 
  VALUE rb_dlamda;
  doublereal *dlamda; 
  VALUE rb_q2;
  doublecomplex *q2; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  doublereal *givnum; 
  VALUE rb_info;
  integer info; 
  VALUE rb_q_out__;
  doublecomplex *q_out__;
  VALUE rb_d_out__;
  doublereal *d_out__;
  integer *indxp;
  integer *indx;

  integer ldq;
  integer n;
  integer ldq2;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, dlamda, q2, w, perm, givptr, givcol, givnum, info, q, d, rho = NumRu::Lapack.zlaed8( qsiz, q, d, rho, cutpnt, z, indxq)\n    or\n  NumRu::Lapack.zlaed8  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z, DLAMDA, Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR, GIVCOL, GIVNUM, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZLAED8 merges the two sets of eigenvalues together into a single\n*  sorted set.  Then it tries to deflate the size of the problem.\n*  There are two ways in which deflation can occur:  when two or more\n*  eigenvalues are close together or if there is a tiny element in the\n*  Z vector.  For each such occurrence the order of the related secular\n*  equation problem is reduced by one.\n*\n\n*  Arguments\n*  =========\n*\n*  K      (output) INTEGER\n*         Contains the number of non-deflated eigenvalues.\n*         This is the order of the related secular equation.\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  QSIZ   (input) INTEGER\n*         The dimension of the unitary matrix used to reduce\n*         the dense or band matrix to tridiagonal form.\n*         QSIZ >= N if ICOMPQ = 1.\n*\n*  Q      (input/output) COMPLEX*16 array, dimension (LDQ,N)\n*         On entry, Q contains the eigenvectors of the partially solved\n*         system which has been previously updated in matrix\n*         multiplies with other partially solved eigensystems.\n*         On exit, Q contains the trailing (N-K) updated eigenvectors\n*         (those which were deflated) in its last N-K columns.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= max( 1, N ).\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension (N)\n*         On entry, D contains the eigenvalues of the two submatrices to\n*         be combined.  On exit, D contains the trailing (N-K) updated\n*         eigenvalues (those which were deflated) sorted into increasing\n*         order.\n*\n*  RHO    (input/output) DOUBLE PRECISION\n*         Contains the off diagonal element associated with the rank-1\n*         cut which originally split the two submatrices which are now\n*         being recombined. RHO is modified during the computation to\n*         the value required by DLAED3.\n*\n*  CUTPNT (input) INTEGER\n*         Contains the location of the last eigenvalue in the leading\n*         sub-matrix.  MIN(1,N) <= CUTPNT <= N.\n*\n*  Z      (input) DOUBLE PRECISION array, dimension (N)\n*         On input this vector contains the updating vector (the last\n*         row of the first sub-eigenvector matrix and the first row of\n*         the second sub-eigenvector matrix).  The contents of Z are\n*         destroyed during the updating process.\n*\n*  DLAMDA (output) DOUBLE PRECISION array, dimension (N)\n*         Contains a copy of the first K eigenvalues which will be used\n*         by DLAED3 to form the secular equation.\n*\n*  Q2     (output) COMPLEX*16 array, dimension (LDQ2,N)\n*         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,\n*         Contains a copy of the first K eigenvectors which will be used\n*         by DLAED7 in a matrix multiply (DGEMM) to update the new\n*         eigenvectors.\n*\n*  LDQ2   (input) INTEGER\n*         The leading dimension of the array Q2.  LDQ2 >= max( 1, N ).\n*\n*  W      (output) DOUBLE PRECISION array, dimension (N)\n*         This will hold the first k values of the final\n*         deflation-altered z-vector and will be passed to DLAED3.\n*\n*  INDXP  (workspace) INTEGER array, dimension (N)\n*         This will contain the permutation used to place deflated\n*         values of D at the end of the array. On output INDXP(1:K)\n*         points to the nondeflated D-values and INDXP(K+1:N)\n*         points to the deflated eigenvalues.\n*\n*  INDX   (workspace) INTEGER array, dimension (N)\n*         This will contain the permutation used to sort the contents of\n*         D into ascending order.\n*\n*  INDXQ  (input) INTEGER array, dimension (N)\n*         This contains the permutation which separately sorts the two\n*         sub-problems in D into ascending order.  Note that elements in\n*         the second half of this permutation must first have CUTPNT\n*         added to their values in order to be accurate.\n*\n*  PERM   (output) INTEGER array, dimension (N)\n*         Contains the permutations (from deflation and sorting) to be\n*         applied to each eigenblock.\n*\n*  GIVPTR (output) INTEGER\n*         Contains the number of Givens rotations which took place in\n*         this subproblem.\n*\n*  GIVCOL (output) INTEGER array, dimension (2, N)\n*         Each pair of numbers indicates a pair of columns to take place\n*         in a Givens rotation.\n*\n*  GIVNUM (output) DOUBLE PRECISION array, dimension (2, N)\n*         Each number indicates the S value to be used in the\n*         corresponding Givens rotation.\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_qsiz = argv[0];
  rb_q = argv[1];
  rb_d = argv[2];
  rb_rho = argv[3];
  rb_cutpnt = argv[4];
  rb_z = argv[5];
  rb_indxq = argv[6];

  qsiz = NUM2INT(rb_qsiz);
  rho = NUM2DBL(rb_rho);
  cutpnt = NUM2INT(rb_cutpnt);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (2th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (2th argument) must be %d", 2);
  ldq = NA_SHAPE0(rb_q);
  n = NA_SHAPE1(rb_q);
  if (NA_TYPE(rb_q) != NA_DCOMPLEX)
    rb_q = na_change_type(rb_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rb_q, doublecomplex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of q");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (6th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 1 of q");
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  if (!NA_IsNArray(rb_indxq))
    rb_raise(rb_eArgError, "indxq (7th argument) must be NArray");
  if (NA_RANK(rb_indxq) != 1)
    rb_raise(rb_eArgError, "rank of indxq (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_indxq) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of indxq must be the same as shape 1 of q");
  if (NA_TYPE(rb_indxq) != NA_LINT)
    rb_indxq = na_change_type(rb_indxq, NA_LINT);
  indxq = NA_PTR_TYPE(rb_indxq, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_dlamda = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dlamda = NA_PTR_TYPE(rb_dlamda, doublereal*);
  ldq2 = MAX( 1, n );
  {
    int shape[2];
    shape[0] = ldq2;
    shape[1] = n;
    rb_q2 = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q2 = NA_PTR_TYPE(rb_q2, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_perm = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  perm = NA_PTR_TYPE(rb_perm, integer*);
  {
    int shape[2];
    shape[0] = 2;
    shape[1] = n;
    rb_givcol = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  {
    int shape[2];
    shape[0] = 2;
    shape[1] = n;
    rb_givnum = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rb_givnum, doublereal*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublecomplex*);
  MEMCPY(q_out__, q, doublecomplex, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  indxp = ALLOC_N(integer, (n));
  indx = ALLOC_N(integer, (n));

  zlaed8_(&k, &n, &qsiz, q, &ldq, d, &rho, &cutpnt, z, dlamda, q2, &ldq2, w, indxp, indx, indxq, perm, &givptr, givcol, givnum, &info);

  free(indxp);
  free(indx);
  rb_k = INT2NUM(k);
  rb_givptr = INT2NUM(givptr);
  rb_info = INT2NUM(info);
  rb_rho = rb_float_new((double)rho);
  return rb_ary_new3(12, rb_k, rb_dlamda, rb_q2, rb_w, rb_perm, rb_givptr, rb_givcol, rb_givnum, rb_info, rb_q, rb_d, rb_rho);
}

void
init_lapack_zlaed8(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaed8", rb_zlaed8, -1);
}
