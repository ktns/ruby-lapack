#include "rb_lapack.h"

extern VOID dsposv_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *work, real *swork, integer *iter, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dsposv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_iter;
  integer iter; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  doublereal *work;
  real *swork;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, iter, info, a = NumRu::Lapack.dsposv( uplo, a, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSPOSV( UPLO, N, NRHS, A, LDA, B, LDB, X, LDX, WORK, SWORK, ITER, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSPOSV computes the solution to a real system of linear equations\n*     A * X = B,\n*  where A is an N-by-N symmetric positive definite matrix and X and B\n*  are N-by-NRHS matrices.\n*\n*  DSPOSV first attempts to factorize the matrix in SINGLE PRECISION\n*  and use this factorization within an iterative refinement procedure\n*  to produce a solution with DOUBLE PRECISION normwise backward error\n*  quality (see below). If the approach fails the method switches to a\n*  DOUBLE PRECISION factorization and solve.\n*\n*  The iterative refinement is not going to be a winning strategy if\n*  the ratio SINGLE PRECISION performance over DOUBLE PRECISION\n*  performance is too small. A reasonable strategy should take the\n*  number of right-hand sides and the size of the matrix into account.\n*  This might be done with a call to ILAENV in the future. Up to now, we\n*  always try iterative refinement.\n*\n*  The iterative refinement process is stopped if\n*      ITER > ITERMAX\n*  or for all the RHS we have:\n*      RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX\n*  where\n*      o ITER is the number of the current iteration in the iterative\n*        refinement process\n*      o RNRM is the infinity-norm of the residual\n*      o XNRM is the infinity-norm of the solution\n*      o ANRM is the infinity-operator-norm of the matrix A\n*      o EPS is the machine epsilon returned by DLAMCH('Epsilon')\n*  The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00\n*  respectively.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array,\n*          dimension (LDA,N)\n*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading\n*          N-by-N upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading N-by-N lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*          On exit, if iterative refinement has been successfully used\n*          (INFO.EQ.0 and ITER.GE.0, see description below), then A is\n*          unchanged, if double precision factorization has been used\n*          (INFO.EQ.0 and ITER.LT.0, see description below), then the\n*          array A contains the factor U or L from the Cholesky\n*          factorization A = U**T*U or A = L*L**T.\n*\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          The N-by-NRHS right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)\n*          If INFO = 0, the N-by-NRHS solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (N,NRHS)\n*          This array is used to hold the residual vectors.\n*\n*  SWORK   (workspace) REAL array, dimension (N*(N+NRHS))\n*          This array is used to use the single precision matrix and the\n*          right-hand sides or solutions in single precision.\n*\n*  ITER    (output) INTEGER\n*          < 0: iterative refinement has failed, double precision\n*               factorization has been performed\n*               -1 : the routine fell back to full precision for\n*                    implementation- or machine-specific reasons\n*               -2 : narrowing the precision induced an overflow,\n*                    the routine fell back to full precision\n*               -3 : failure of SPOTRF\n*               -31: stop the iterative refinement after the 30th\n*                    iterations\n*          > 0: iterative refinement has been sucessfully used.\n*               Returns the number of iterations\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the leading minor of order i of (DOUBLE\n*                PRECISION) A is not positive definite, so the\n*                factorization could not be completed, and the solution\n*                has not been computed.\n*\n*  =========\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, iter, info, a = NumRu::Lapack.dsposv( uplo, a, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_uplo = argv[0];
  rblapack_a = argv[1];
  rblapack_b = argv[2];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rblapack_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  work = ALLOC_N(doublereal, (n)*(nrhs));
  swork = ALLOC_N(real, (n*(n+nrhs)));

  dsposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, work, swork, &iter, &info);

  free(work);
  free(swork);
  rblapack_iter = INT2NUM(iter);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_x, rblapack_iter, rblapack_info, rblapack_a);
}

void
init_lapack_dsposv(VALUE mLapack){
  rb_define_module_function(mLapack, "dsposv", rblapack_dsposv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
