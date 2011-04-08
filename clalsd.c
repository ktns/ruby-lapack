#include "rb_lapack.h"

extern VOID clalsd_(char *uplo, integer *smlsiz, integer *n, integer *nrhs, real *d, real *e, complex *b, integer *ldb, real *rcond, integer *rank, complex *work, real *rwork, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clalsd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_smlsiz;
  integer smlsiz; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_rcond;
  real rcond; 
  VALUE rblapack_rank;
  integer rank; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_e_out__;
  real *e_out__;
  VALUE rblapack_b_out__;
  complex *b_out__;
  integer nlvl;
  complex *work;
  real *rwork;
  integer *iwork;

  integer n;
  integer ldb;
  integer nrhs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rank, info, d, e, b = NumRu::Lapack.clalsd( uplo, smlsiz, d, e, b, rcond, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, RANK, WORK, RWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CLALSD uses the singular value decomposition of A to solve the least\n*  squares problem of finding X to minimize the Euclidean norm of each\n*  column of A*X-B, where A is N-by-N upper bidiagonal, and X and B\n*  are N-by-NRHS. The solution X overwrites B.\n*\n*  The singular values of A smaller than RCOND times the largest\n*  singular value are treated as zero in solving the least squares\n*  problem; in this case a minimum norm solution is returned.\n*  The actual singular values are returned in D in ascending order.\n*\n*  This code makes very mild assumptions about floating point\n*  arithmetic. It will work on machines with a guard digit in\n*  add/subtract, or on those binary machines without guard digits\n*  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.\n*  It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO   (input) CHARACTER*1\n*         = 'U': D and E define an upper bidiagonal matrix.\n*         = 'L': D and E define a  lower bidiagonal matrix.\n*\n*  SMLSIZ (input) INTEGER\n*         The maximum size of the subproblems at the bottom of the\n*         computation tree.\n*\n*  N      (input) INTEGER\n*         The dimension of the  bidiagonal matrix.  N >= 0.\n*\n*  NRHS   (input) INTEGER\n*         The number of columns of B. NRHS must be at least 1.\n*\n*  D      (input/output) REAL array, dimension (N)\n*         On entry D contains the main diagonal of the bidiagonal\n*         matrix. On exit, if INFO = 0, D contains its singular values.\n*\n*  E      (input/output) REAL array, dimension (N-1)\n*         Contains the super-diagonal entries of the bidiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  B      (input/output) COMPLEX array, dimension (LDB,NRHS)\n*         On input, B contains the right hand sides of the least\n*         squares problem. On output, B contains the solution X.\n*\n*  LDB    (input) INTEGER\n*         The leading dimension of B in the calling subprogram.\n*         LDB must be at least max(1,N).\n*\n*  RCOND  (input) REAL\n*         The singular values of A less than or equal to RCOND times\n*         the largest singular value are treated as zero in solving\n*         the least squares problem. If RCOND is negative,\n*         machine precision is used instead.\n*         For example, if diag(S)*X=B were the least squares problem,\n*         where diag(S) is a diagonal matrix of singular values, the\n*         solution would be X(i) = B(i) / S(i) if S(i) is greater than\n*         RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to\n*         RCOND*max(S).\n*\n*  RANK   (output) INTEGER\n*         The number of singular values of A greater than RCOND times\n*         the largest singular value.\n*\n*  WORK   (workspace) COMPLEX array, dimension (N * NRHS).\n*\n*  RWORK  (workspace) REAL array, dimension at least\n*         (9*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS +\n*         MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS ),\n*         where\n*         NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )\n*\n*  IWORK  (workspace) INTEGER array, dimension (3*N*NLVL + 11*N).\n*\n*  INFO   (output) INTEGER\n*         = 0:  successful exit.\n*         < 0:  if INFO = -i, the i-th argument had an illegal value.\n*         > 0:  The algorithm failed to compute a singular value while\n*               working on the submatrix lying in rows and columns\n*               INFO/(N+1) through MOD(INFO,N+1).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Ren-Cang Li, Computer Science Division, University of\n*       California at Berkeley, USA\n*     Osni Marques, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rank, info, d, e, b = NumRu::Lapack.clalsd( uplo, smlsiz, d, e, b, rcond, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_uplo = argv[0];
  rblapack_smlsiz = argv[1];
  rblapack_d = argv[2];
  rblapack_e = argv[3];
  rblapack_b = argv[4];
  rblapack_rcond = argv[5];
  if (rb_options != Qnil) {
  }

  rcond = (real)NUM2DBL(rblapack_rcond);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  smlsiz = NUM2INT(rblapack_smlsiz);
  uplo = StringValueCStr(rblapack_uplo)[0];
  nlvl = ( (int)( log(((double)n)/(smlsiz+1))/log(2.0) ) ) + 1;
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  work = ALLOC_N(complex, (n * nrhs));
  rwork = ALLOC_N(real, (9*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1)));
  iwork = ALLOC_N(integer, (3*n*nlvl + 11*n));

  clalsd_(&uplo, &smlsiz, &n, &nrhs, d, e, b, &ldb, &rcond, &rank, work, rwork, iwork, &info);

  free(work);
  free(rwork);
  free(iwork);
  rblapack_rank = INT2NUM(rank);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_rank, rblapack_info, rblapack_d, rblapack_e, rblapack_b);
}

void
init_lapack_clalsd(VALUE mLapack){
  rb_define_module_function(mLapack, "clalsd", rblapack_clalsd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
