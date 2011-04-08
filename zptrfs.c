#include "rb_lapack.h"

extern VOID zptrfs_(char *uplo, integer *n, integer *nrhs, doublereal *d, doublecomplex *e, doublereal *df, doublecomplex *ef, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zptrfs(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublecomplex *e; 
  VALUE rblapack_df;
  doublereal *df; 
  VALUE rblapack_ef;
  doublecomplex *ef; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_x;
  doublecomplex *x; 
  VALUE rblapack_ferr;
  doublereal *ferr; 
  VALUE rblapack_berr;
  doublereal *berr; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_x_out__;
  doublecomplex *x_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.zptrfs( uplo, d, e, df, ef, b, x, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZPTRFS improves the computed solution to a system of linear\n*  equations when the coefficient matrix is Hermitian positive definite\n*  and tridiagonal, and provides error bounds and backward error\n*  estimates for the solution.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the superdiagonal or the subdiagonal of the\n*          tridiagonal matrix A is stored and the form of the\n*          factorization:\n*          = 'U':  E is the superdiagonal of A, and A = U**H*D*U;\n*          = 'L':  E is the subdiagonal of A, and A = L*D*L**H.\n*          (The two forms are equivalent if A is real.)\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The n real diagonal elements of the tridiagonal matrix A.\n*\n*  E       (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) off-diagonal elements of the tridiagonal matrix A\n*          (see UPLO).\n*\n*  DF      (input) DOUBLE PRECISION array, dimension (N)\n*          The n diagonal elements of the diagonal matrix D from\n*          the factorization computed by ZPTTRF.\n*\n*  EF      (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) off-diagonal elements of the unit bidiagonal\n*          factor U or L from the factorization computed by ZPTTRF\n*          (see UPLO).\n*\n*  B       (input) COMPLEX*16 array, dimension (LDB,NRHS)\n*          The right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (input/output) COMPLEX*16 array, dimension (LDX,NRHS)\n*          On entry, the solution matrix X, as computed by ZPTTRS.\n*          On exit, the improved solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).\n*\n*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n*  Internal Parameters\n*  ===================\n*\n*  ITMAX is the maximum number of steps of iterative refinement.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.zptrfs( uplo, d, e, df, ef, b, x, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_uplo = argv[0];
  rblapack_d = argv[1];
  rblapack_e = argv[2];
  rblapack_df = argv[3];
  rblapack_ef = argv[4];
  rblapack_b = argv[5];
  rblapack_x = argv[6];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (7th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 2)
    rb_raise(rb_eArgError, "rank of x (7th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_x);
  ldx = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_DCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  if (!NA_IsNArray(rblapack_df))
    rb_raise(rb_eArgError, "df (4th argument) must be NArray");
  if (NA_RANK(rblapack_df) != 1)
    rb_raise(rb_eArgError, "rank of df (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_df) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of df must be the same as shape 0 of d");
  if (NA_TYPE(rblapack_df) != NA_DFLOAT)
    rblapack_df = na_change_type(rblapack_df, NA_DFLOAT);
  df = NA_PTR_TYPE(rblapack_df, doublereal*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_DCOMPLEX)
    rblapack_e = na_change_type(rblapack_e, NA_DCOMPLEX);
  e = NA_PTR_TYPE(rblapack_e, doublecomplex*);
  if (!NA_IsNArray(rblapack_ef))
    rb_raise(rb_eArgError, "ef (5th argument) must be NArray");
  if (NA_RANK(rblapack_ef) != 1)
    rb_raise(rb_eArgError, "rank of ef (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ef) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ef must be %d", n-1);
  if (NA_TYPE(rblapack_ef) != NA_DCOMPLEX)
    rblapack_ef = na_change_type(rblapack_ef, NA_DCOMPLEX);
  ef = NA_PTR_TYPE(rblapack_ef, doublecomplex*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_ferr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rblapack_ferr, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rblapack_berr, doublereal*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rblapack_x_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  work = ALLOC_N(doublecomplex, (n));
  rwork = ALLOC_N(doublereal, (n));

  zptrfs_(&uplo, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_ferr, rblapack_berr, rblapack_info, rblapack_x);
}

void
init_lapack_zptrfs(VALUE mLapack){
  rb_define_module_function(mLapack, "zptrfs", rblapack_zptrfs, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
