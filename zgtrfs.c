#include "rb_lapack.h"

extern VOID zgtrfs_(char *trans, integer *n, integer *nrhs, doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgtrfs(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_dl;
  doublecomplex *dl; 
  VALUE rblapack_d;
  doublecomplex *d; 
  VALUE rblapack_du;
  doublecomplex *du; 
  VALUE rblapack_dlf;
  doublecomplex *dlf; 
  VALUE rblapack_df;
  doublecomplex *df; 
  VALUE rblapack_duf;
  doublecomplex *duf; 
  VALUE rblapack_du2;
  doublecomplex *du2; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
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
      printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.zgtrfs( trans, dl, d, du, dlf, df, duf, du2, ipiv, b, x, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGTRFS improves the computed solution to a system of linear\n*  equations when the coefficient matrix is tridiagonal, and provides\n*  error bounds and backward error estimates for the solution.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose)\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  DL      (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) subdiagonal elements of A.\n*\n*  D       (input) COMPLEX*16 array, dimension (N)\n*          The diagonal elements of A.\n*\n*  DU      (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) superdiagonal elements of A.\n*\n*  DLF     (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) multipliers that define the matrix L from the\n*          LU factorization of A as computed by ZGTTRF.\n*\n*  DF      (input) COMPLEX*16 array, dimension (N)\n*          The n diagonal elements of the upper triangular matrix U from\n*          the LU factorization of A.\n*\n*  DUF     (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) elements of the first superdiagonal of U.\n*\n*  DU2     (input) COMPLEX*16 array, dimension (N-2)\n*          The (n-2) elements of the second superdiagonal of U.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= n, row i of the matrix was\n*          interchanged with row IPIV(i).  IPIV(i) will always be either\n*          i or i+1; IPIV(i) = i indicates a row interchange was not\n*          required.\n*\n*  B       (input) COMPLEX*16 array, dimension (LDB,NRHS)\n*          The right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (input/output) COMPLEX*16 array, dimension (LDX,NRHS)\n*          On entry, the solution matrix X, as computed by ZGTTRS.\n*          On exit, the improved solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (2*N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n*  Internal Parameters\n*  ===================\n*\n*  ITMAX is the maximum number of steps of iterative refinement.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.zgtrfs( trans, dl, d, du, dlf, df, duf, du2, ipiv, b, x, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_x) != NA_DCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (10th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (10th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  if (!NA_IsNArray(rblapack_df))
    rb_raise(rb_eArgError, "df (6th argument) must be NArray");
  if (NA_RANK(rblapack_df) != 1)
    rb_raise(rb_eArgError, "rank of df (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_df) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of df must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_df) != NA_DCOMPLEX)
    rblapack_df = na_change_type(rblapack_df, NA_DCOMPLEX);
  df = NA_PTR_TYPE(rblapack_df, doublecomplex*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_d) != NA_DCOMPLEX)
    rblapack_d = na_change_type(rblapack_d, NA_DCOMPLEX);
  d = NA_PTR_TYPE(rblapack_d, doublecomplex*);
  if (!NA_IsNArray(rblapack_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rblapack_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rblapack_du) != NA_DCOMPLEX)
    rblapack_du = na_change_type(rblapack_du, NA_DCOMPLEX);
  du = NA_PTR_TYPE(rblapack_du, doublecomplex*);
  if (!NA_IsNArray(rblapack_dlf))
    rb_raise(rb_eArgError, "dlf (5th argument) must be NArray");
  if (NA_RANK(rblapack_dlf) != 1)
    rb_raise(rb_eArgError, "rank of dlf (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dlf) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dlf must be %d", n-1);
  if (NA_TYPE(rblapack_dlf) != NA_DCOMPLEX)
    rblapack_dlf = na_change_type(rblapack_dlf, NA_DCOMPLEX);
  dlf = NA_PTR_TYPE(rblapack_dlf, doublecomplex*);
  if (!NA_IsNArray(rblapack_dl))
    rb_raise(rb_eArgError, "dl (2th argument) must be NArray");
  if (NA_RANK(rblapack_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rblapack_dl) != NA_DCOMPLEX)
    rblapack_dl = na_change_type(rblapack_dl, NA_DCOMPLEX);
  dl = NA_PTR_TYPE(rblapack_dl, doublecomplex*);
  if (!NA_IsNArray(rblapack_duf))
    rb_raise(rb_eArgError, "duf (7th argument) must be NArray");
  if (NA_RANK(rblapack_duf) != 1)
    rb_raise(rb_eArgError, "rank of duf (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_duf) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of duf must be %d", n-1);
  if (NA_TYPE(rblapack_duf) != NA_DCOMPLEX)
    rblapack_duf = na_change_type(rblapack_duf, NA_DCOMPLEX);
  duf = NA_PTR_TYPE(rblapack_duf, doublecomplex*);
  if (!NA_IsNArray(rblapack_du2))
    rb_raise(rb_eArgError, "du2 (8th argument) must be NArray");
  if (NA_RANK(rblapack_du2) != 1)
    rb_raise(rb_eArgError, "rank of du2 (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du2) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of du2 must be %d", n-2);
  if (NA_TYPE(rblapack_du2) != NA_DCOMPLEX)
    rblapack_du2 = na_change_type(rblapack_du2, NA_DCOMPLEX);
  du2 = NA_PTR_TYPE(rblapack_du2, doublecomplex*);
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
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (n));

  zgtrfs_(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_ferr, rblapack_berr, rblapack_info, rblapack_x);
}

void
init_lapack_zgtrfs(VALUE mLapack){
  rb_define_module_function(mLapack, "zgtrfs", rblapack_zgtrfs, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
