#include "rb_lapack.h"

static VALUE
rb_zherfs(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_af;
  doublecomplex *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_ferr;
  doublereal *ferr; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x_out__;
  doublecomplex *x_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.zherfs( uplo, a, af, ipiv, b, x)\n    or\n  NumRu::Lapack.zherfs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHERFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZHERFS improves the computed solution to a system of linear\n*  equations when the coefficient matrix is Hermitian indefinite, and\n*  provides error bounds and backward error estimates for the solution.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA,N)\n*          The Hermitian matrix A.  If UPLO = 'U', the leading N-by-N\n*          upper triangular part of A contains the upper triangular part\n*          of the matrix A, and the strictly lower triangular part of A\n*          is not referenced.  If UPLO = 'L', the leading N-by-N lower\n*          triangular part of A contains the lower triangular part of\n*          the matrix A, and the strictly upper triangular part of A is\n*          not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  AF      (input) COMPLEX*16 array, dimension (LDAF,N)\n*          The factored form of the matrix A.  AF contains the block\n*          diagonal matrix D and the multipliers used to obtain the\n*          factor U or L from the factorization A = U*D*U**H or\n*          A = L*D*L**H as computed by ZHETRF.\n*\n*  LDAF    (input) INTEGER\n*          The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D\n*          as determined by ZHETRF.\n*\n*  B       (input) COMPLEX*16 array, dimension (LDB,NRHS)\n*          The right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (input/output) COMPLEX*16 array, dimension (LDX,NRHS)\n*          On entry, the solution matrix X, as computed by ZHETRS.\n*          On exit, the improved solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (2*N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n*  Internal Parameters\n*  ===================\n*\n*  ITMAX is the maximum number of steps of iterative refinement.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];
  rb_ipiv = argv[3];
  rb_b = argv[4];
  rb_x = argv[5];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  ldaf = NA_SHAPE0(rb_af);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  if (NA_TYPE(rb_af) != NA_DCOMPLEX)
    rb_af = na_change_type(rb_af, NA_DCOMPLEX);
  af = NA_PTR_TYPE(rb_af, doublecomplex*);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ipiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be the same as shape 1 of a");
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 2);
  ldx = NA_SHAPE0(rb_x);
  if (NA_SHAPE1(rb_x) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of x must be the same as shape 1 of b");
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_ferr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rb_ferr, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, doublereal*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (n));

  zherfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_ferr, rb_berr, rb_info, rb_x);
}

void
init_lapack_zherfs(VALUE mLapack){
  rb_define_module_function(mLapack, "zherfs", rb_zherfs, -1);
}
