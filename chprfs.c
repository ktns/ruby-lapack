#include "rb_lapack.h"

extern VOID chprfs_(char *uplo, integer *n, integer *nrhs, complex *ap, complex *afp, integer *ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_chprfs(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ap;
  complex *ap; 
  VALUE rblapack_afp;
  complex *afp; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_x;
  complex *x; 
  VALUE rblapack_ferr;
  real *ferr; 
  VALUE rblapack_berr;
  real *berr; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_x_out__;
  complex *x_out__;
  complex *work;
  real *rwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.chprfs( uplo, ap, afp, ipiv, b, x, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CHPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CHPRFS improves the computed solution to a system of linear\n*  equations when the coefficient matrix is Hermitian indefinite\n*  and packed, and provides error bounds and backward error estimates\n*  for the solution.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  AP      (input) COMPLEX array, dimension (N*(N+1)/2)\n*          The upper or lower triangle of the Hermitian matrix A, packed\n*          columnwise in a linear array.  The j-th column of A is stored\n*          in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*\n*  AFP     (input) COMPLEX array, dimension (N*(N+1)/2)\n*          The factored form of the matrix A.  AFP contains the block\n*          diagonal matrix D and the multipliers used to obtain the\n*          factor U or L from the factorization A = U*D*U**H or\n*          A = L*D*L**H as computed by CHPTRF, stored as a packed\n*          triangular matrix.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D\n*          as determined by CHPTRF.\n*\n*  B       (input) COMPLEX array, dimension (LDB,NRHS)\n*          The right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (input/output) COMPLEX array, dimension (LDX,NRHS)\n*          On entry, the solution matrix X, as computed by CHPTRS.\n*          On exit, the improved solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  FERR    (output) REAL array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) REAL array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  RWORK   (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n*  Internal Parameters\n*  ===================\n*\n*  ITMAX is the maximum number of steps of iterative refinement.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.chprfs( uplo, ap, afp, ipiv, b, x, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_uplo = argv[0];
  rblapack_ap = argv[1];
  rblapack_afp = argv[2];
  rblapack_ipiv = argv[3];
  rblapack_b = argv[4];
  rblapack_x = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 2)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_x);
  ldx = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_SCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, complex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_ap) != NA_SCOMPLEX)
    rblapack_ap = na_change_type(rblapack_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rblapack_ap, complex*);
  if (!NA_IsNArray(rblapack_afp))
    rb_raise(rb_eArgError, "afp (3th argument) must be NArray");
  if (NA_RANK(rblapack_afp) != 1)
    rb_raise(rb_eArgError, "rank of afp (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_afp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of afp must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_afp) != NA_SCOMPLEX)
    rblapack_afp = na_change_type(rblapack_afp, NA_SCOMPLEX);
  afp = NA_PTR_TYPE(rblapack_afp, complex*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_ferr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rblapack_ferr, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rblapack_berr, real*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rblapack_x_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, complex*);
  MEMCPY(x_out__, x, complex, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (n));

  chprfs_(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_ferr, rblapack_berr, rblapack_info, rblapack_x);
}

void
init_lapack_chprfs(VALUE mLapack){
  rb_define_module_function(mLapack, "chprfs", rblapack_chprfs, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
