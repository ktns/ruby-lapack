#include "rb_lapack.h"

static VALUE
rb_sptrfs(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_df;
  real *df; 
  VALUE rb_ef;
  real *ef; 
  VALUE rb_b;
  real *b; 
  VALUE rb_x;
  real *x; 
  VALUE rb_ferr;
  real *ferr; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x_out__;
  real *x_out__;
  real *work;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.sptrfs( d, e, df, ef, b, x)\n    or\n  NumRu::Lapack.sptrfs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPTRFS improves the computed solution to a system of linear\n*  equations when the coefficient matrix is symmetric positive definite\n*  and tridiagonal, and provides error bounds and backward error\n*  estimates for the solution.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  D       (input) REAL array, dimension (N)\n*          The n diagonal elements of the tridiagonal matrix A.\n*\n*  E       (input) REAL array, dimension (N-1)\n*          The (n-1) subdiagonal elements of the tridiagonal matrix A.\n*\n*  DF      (input) REAL array, dimension (N)\n*          The n diagonal elements of the diagonal matrix D from the\n*          factorization computed by SPTTRF.\n*\n*  EF      (input) REAL array, dimension (N-1)\n*          The (n-1) subdiagonal elements of the unit bidiagonal factor\n*          L from the factorization computed by SPTTRF.\n*\n*  B       (input) REAL array, dimension (LDB,NRHS)\n*          The right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (input/output) REAL array, dimension (LDX,NRHS)\n*          On entry, the solution matrix X, as computed by SPTTRS.\n*          On exit, the improved solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  FERR    (output) REAL array, dimension (NRHS)\n*          The forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).\n*\n*  BERR    (output) REAL array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) REAL array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n*  Internal Parameters\n*  ===================\n*\n*  ITMAX is the maximum number of steps of iterative refinement.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_d = argv[0];
  rb_e = argv[1];
  rb_df = argv[2];
  rb_ef = argv[3];
  rb_b = argv[4];
  rb_x = argv[5];

  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  if (!NA_IsNArray(rb_df))
    rb_raise(rb_eArgError, "df (3th argument) must be NArray");
  if (NA_RANK(rb_df) != 1)
    rb_raise(rb_eArgError, "rank of df (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_df) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of df must be the same as shape 0 of d");
  if (NA_TYPE(rb_df) != NA_SFLOAT)
    rb_df = na_change_type(rb_df, NA_SFLOAT);
  df = NA_PTR_TYPE(rb_df, real*);
  if (!NA_IsNArray(rb_ef))
    rb_raise(rb_eArgError, "ef (4th argument) must be NArray");
  if (NA_RANK(rb_ef) != 1)
    rb_raise(rb_eArgError, "rank of ef (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ef) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ef must be %d", n-1);
  if (NA_TYPE(rb_ef) != NA_SFLOAT)
    rb_ef = na_change_type(rb_ef, NA_SFLOAT);
  ef = NA_PTR_TYPE(rb_ef, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 2);
  ldx = NA_SHAPE0(rb_x);
  if (NA_SHAPE1(rb_x) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of x must be the same as shape 1 of b");
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_ferr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rb_ferr, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, real*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  work = ALLOC_N(real, (2*n));

  sptrfs_(&n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_ferr, rb_berr, rb_info, rb_x);
}

void
init_lapack_sptrfs(VALUE mLapack){
  rb_define_module_function(mLapack, "sptrfs", rb_sptrfs, -1);
}
