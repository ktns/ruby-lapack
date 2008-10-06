#include "rb_lapack.h"

static VALUE
rb_slalsd(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_smlsiz;
  integer smlsiz; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_b;
  real *b; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  VALUE rb_b_out__;
  real *b_out__;
  real *work;
  integer *iwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer nlvl;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rank, info, d, e, b = NumRu::Lapack.slalsd( uplo, smlsiz, d, e, b, rcond)\n    or\n  NumRu::Lapack.slalsd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, RANK, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLALSD uses the singular value decomposition of A to solve the least\n*  squares problem of finding X to minimize the Euclidean norm of each\n*  column of A*X-B, where A is N-by-N upper bidiagonal, and X and B\n*  are N-by-NRHS. The solution X overwrites B.\n*\n*  The singular values of A smaller than RCOND times the largest\n*  singular value are treated as zero in solving the least squares\n*  problem; in this case a minimum norm solution is returned.\n*  The actual singular values are returned in D in ascending order.\n*\n*  This code makes very mild assumptions about floating point\n*  arithmetic. It will work on machines with a guard digit in\n*  add/subtract, or on those binary machines without guard digits\n*  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.\n*  It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO   (input) CHARACTER*1\n*         = 'U': D and E define an upper bidiagonal matrix.\n*         = 'L': D and E define a  lower bidiagonal matrix.\n*\n*  SMLSIZ (input) INTEGER\n*         The maximum size of the subproblems at the bottom of the\n*         computation tree.\n*\n*  N      (input) INTEGER\n*         The dimension of the  bidiagonal matrix.  N >= 0.\n*\n*  NRHS   (input) INTEGER\n*         The number of columns of B. NRHS must be at least 1.\n*\n*  D      (input/output) REAL array, dimension (N)\n*         On entry D contains the main diagonal of the bidiagonal\n*         matrix. On exit, if INFO = 0, D contains its singular values.\n*\n*  E      (input/output) REAL array, dimension (N-1)\n*         Contains the super-diagonal entries of the bidiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  B      (input/output) REAL array, dimension (LDB,NRHS)\n*         On input, B contains the right hand sides of the least\n*         squares problem. On output, B contains the solution X.\n*\n*  LDB    (input) INTEGER\n*         The leading dimension of B in the calling subprogram.\n*         LDB must be at least max(1,N).\n*\n*  RCOND  (input) REAL\n*         The singular values of A less than or equal to RCOND times\n*         the largest singular value are treated as zero in solving\n*         the least squares problem. If RCOND is negative,\n*         machine precision is used instead.\n*         For example, if diag(S)*X=B were the least squares problem,\n*         where diag(S) is a diagonal matrix of singular values, the\n*         solution would be X(i) = B(i) / S(i) if S(i) is greater than\n*         RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to\n*         RCOND*max(S).\n*\n*  RANK   (output) INTEGER\n*         The number of singular values of A greater than RCOND times\n*         the largest singular value.\n*\n*  WORK   (workspace) REAL array, dimension at least\n*         (9*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2),\n*         where NLVL = max(0, INT(log_2 (N/(SMLSIZ+1))) + 1).\n*\n*  IWORK  (workspace) INTEGER array, dimension at least\n*         (3*N*NLVL + 11*N)\n*\n*  INFO   (output) INTEGER\n*         = 0:  successful exit.\n*         < 0:  if INFO = -i, the i-th argument had an illegal value.\n*         > 0:  The algorithm failed to compute an singular value while\n*               working on the submatrix lying in rows and columns\n*               INFO/(N+1) through MOD(INFO,N+1).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Ren-Cang Li, Computer Science Division, University of\n*       California at Berkeley, USA\n*     Osni Marques, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_smlsiz = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_b = argv[4];
  rb_rcond = argv[5];

  uplo = StringValueCStr(rb_uplo)[0];
  smlsiz = NUM2INT(rb_smlsiz);
  rcond = (real)NUM2DBL(rb_rcond);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
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
    int shape[1];
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  nlvl = MAX(0, (int)(1.0/log(2.0)*log((double)n/(smlsiz+1))) + 1);
  work = ALLOC_N(real, (9*n + 2*n*smlsiz + 8*n*nlvl + n*nrhs + pow(smlsiz+1,2)));
  iwork = ALLOC_N(integer, (3*n*nlvl + 11*n));

  slalsd_(&uplo, &smlsiz, &n, &nrhs, d, e, b, &ldb, &rcond, &rank, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_rank, rb_info, rb_d, rb_e, rb_b);
}

void
init_lapack_slalsd(VALUE mLapack){
  rb_define_module_function(mLapack, "slalsd", rb_slalsd, -1);
}
