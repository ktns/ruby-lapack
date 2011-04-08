#include "rb_lapack.h"

extern VOID cgtrfs_(char *trans, integer *n, integer *nrhs, complex *dl, complex *d, complex *du, complex *dlf, complex *df, complex *duf, complex *du2, integer *ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgtrfs(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_dl;
  complex *dl; 
  VALUE rblapack_d;
  complex *d; 
  VALUE rblapack_du;
  complex *du; 
  VALUE rblapack_dlf;
  complex *dlf; 
  VALUE rblapack_df;
  complex *df; 
  VALUE rblapack_duf;
  complex *duf; 
  VALUE rblapack_du2;
  complex *du2; 
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
      printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.cgtrfs( trans, dl, d, du, dlf, df, duf, du2, ipiv, b, x, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGTRFS improves the computed solution to a system of linear\n*  equations when the coefficient matrix is tridiagonal, and provides\n*  error bounds and backward error estimates for the solution.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose)\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  DL      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) subdiagonal elements of A.\n*\n*  D       (input) COMPLEX array, dimension (N)\n*          The diagonal elements of A.\n*\n*  DU      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) superdiagonal elements of A.\n*\n*  DLF     (input) COMPLEX array, dimension (N-1)\n*          The (n-1) multipliers that define the matrix L from the\n*          LU factorization of A as computed by CGTTRF.\n*\n*  DF      (input) COMPLEX array, dimension (N)\n*          The n diagonal elements of the upper triangular matrix U from\n*          the LU factorization of A.\n*\n*  DUF     (input) COMPLEX array, dimension (N-1)\n*          The (n-1) elements of the first superdiagonal of U.\n*\n*  DU2     (input) COMPLEX array, dimension (N-2)\n*          The (n-2) elements of the second superdiagonal of U.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= n, row i of the matrix was\n*          interchanged with row IPIV(i).  IPIV(i) will always be either\n*          i or i+1; IPIV(i) = i indicates a row interchange was not\n*          required.\n*\n*  B       (input) COMPLEX array, dimension (LDB,NRHS)\n*          The right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (input/output) COMPLEX array, dimension (LDX,NRHS)\n*          On entry, the solution matrix X, as computed by CGTTRS.\n*          On exit, the improved solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  FERR    (output) REAL array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) REAL array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  RWORK   (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n*  Internal Parameters\n*  ===================\n*\n*  ITMAX is the maximum number of steps of iterative refinement.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.cgtrfs( trans, dl, d, du, dlf, df, duf, du2, ipiv, b, x, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rblapack_trans = argv[0];
  rblapack_dl = argv[1];
  rblapack_d = argv[2];
  rblapack_du = argv[3];
  rblapack_dlf = argv[4];
  rblapack_df = argv[5];
  rblapack_duf = argv[6];
  rblapack_du2 = argv[7];
  rblapack_ipiv = argv[8];
  rblapack_b = argv[9];
  rblapack_x = argv[10];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (9th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (9th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  trans = StringValueCStr(rblapack_trans)[0];
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (11th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 2)
    rb_raise(rb_eArgError, "rank of x (11th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_x);
  ldx = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_SCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, complex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (10th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (10th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  if (!NA_IsNArray(rblapack_df))
    rb_raise(rb_eArgError, "df (6th argument) must be NArray");
  if (NA_RANK(rblapack_df) != 1)
    rb_raise(rb_eArgError, "rank of df (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_df) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of df must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_df) != NA_SCOMPLEX)
    rblapack_df = na_change_type(rblapack_df, NA_SCOMPLEX);
  df = NA_PTR_TYPE(rblapack_df, complex*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_d) != NA_SCOMPLEX)
    rblapack_d = na_change_type(rblapack_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rblapack_d, complex*);
  if (!NA_IsNArray(rblapack_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rblapack_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rblapack_du) != NA_SCOMPLEX)
    rblapack_du = na_change_type(rblapack_du, NA_SCOMPLEX);
  du = NA_PTR_TYPE(rblapack_du, complex*);
  if (!NA_IsNArray(rblapack_dlf))
    rb_raise(rb_eArgError, "dlf (5th argument) must be NArray");
  if (NA_RANK(rblapack_dlf) != 1)
    rb_raise(rb_eArgError, "rank of dlf (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dlf) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dlf must be %d", n-1);
  if (NA_TYPE(rblapack_dlf) != NA_SCOMPLEX)
    rblapack_dlf = na_change_type(rblapack_dlf, NA_SCOMPLEX);
  dlf = NA_PTR_TYPE(rblapack_dlf, complex*);
  if (!NA_IsNArray(rblapack_dl))
    rb_raise(rb_eArgError, "dl (2th argument) must be NArray");
  if (NA_RANK(rblapack_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rblapack_dl) != NA_SCOMPLEX)
    rblapack_dl = na_change_type(rblapack_dl, NA_SCOMPLEX);
  dl = NA_PTR_TYPE(rblapack_dl, complex*);
  if (!NA_IsNArray(rblapack_duf))
    rb_raise(rb_eArgError, "duf (7th argument) must be NArray");
  if (NA_RANK(rblapack_duf) != 1)
    rb_raise(rb_eArgError, "rank of duf (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_duf) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of duf must be %d", n-1);
  if (NA_TYPE(rblapack_duf) != NA_SCOMPLEX)
    rblapack_duf = na_change_type(rblapack_duf, NA_SCOMPLEX);
  duf = NA_PTR_TYPE(rblapack_duf, complex*);
  if (!NA_IsNArray(rblapack_du2))
    rb_raise(rb_eArgError, "du2 (8th argument) must be NArray");
  if (NA_RANK(rblapack_du2) != 1)
    rb_raise(rb_eArgError, "rank of du2 (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du2) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of du2 must be %d", n-2);
  if (NA_TYPE(rblapack_du2) != NA_SCOMPLEX)
    rblapack_du2 = na_change_type(rblapack_du2, NA_SCOMPLEX);
  du2 = NA_PTR_TYPE(rblapack_du2, complex*);
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

  cgtrfs_(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_ferr, rblapack_berr, rblapack_info, rblapack_x);
}

void
init_lapack_cgtrfs(VALUE mLapack){
  rb_define_module_function(mLapack, "cgtrfs", rblapack_cgtrfs, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
