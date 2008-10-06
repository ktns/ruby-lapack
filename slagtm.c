#include "rb_lapack.h"

static VALUE
rb_slagtm(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_alpha;
  real alpha; 
  VALUE rb_dl;
  real *dl; 
  VALUE rb_d;
  real *d; 
  VALUE rb_du;
  real *du; 
  VALUE rb_x;
  real *x; 
  VALUE rb_beta;
  real beta; 
  VALUE rb_b;
  real *b; 
  VALUE rb_b_out__;
  real *b_out__;

  integer n;
  integer ldx;
  integer nrhs;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  b = NumRu::Lapack.slagtm( trans, alpha, dl, d, du, x, beta, b)\n    or\n  NumRu::Lapack.slagtm  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, B, LDB )\n\n*  Purpose\n*  =======\n*\n*  SLAGTM performs a matrix-vector product of the form\n*\n*     B := alpha * A * X + beta * B\n*\n*  where A is a tridiagonal matrix of order N, B and X are N by NRHS\n*  matrices, and alpha and beta are real scalars, each of which may be\n*  0., 1., or -1.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the operation applied to A.\n*          = 'N':  No transpose, B := alpha * A * X + beta * B\n*          = 'T':  Transpose,    B := alpha * A'* X + beta * B\n*          = 'C':  Conjugate transpose = Transpose\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices X and B.\n*\n*  ALPHA   (input) REAL\n*          The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise,\n*          it is assumed to be 0.\n*\n*  DL      (input) REAL array, dimension (N-1)\n*          The (n-1) sub-diagonal elements of T.\n*\n*  D       (input) REAL array, dimension (N)\n*          The diagonal elements of T.\n*\n*  DU      (input) REAL array, dimension (N-1)\n*          The (n-1) super-diagonal elements of T.\n*\n*  X       (input) REAL array, dimension (LDX,NRHS)\n*          The N by NRHS matrix X.\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(N,1).\n*\n*  BETA    (input) REAL\n*          The scalar beta.  BETA must be 0., 1., or -1.; otherwise,\n*          it is assumed to be 1.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the N by NRHS matrix B.\n*          On exit, B is overwritten by the matrix expression\n*          B := alpha * A * X + beta * B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(N,1).\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_trans = argv[0];
  rb_alpha = argv[1];
  rb_dl = argv[2];
  rb_d = argv[3];
  rb_du = argv[4];
  rb_x = argv[5];
  rb_beta = argv[6];
  rb_b = argv[7];

  trans = StringValueCStr(rb_trans)[0];
  alpha = (real)NUM2DBL(rb_alpha);
  beta = (real)NUM2DBL(rb_beta);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (5th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SFLOAT)
    rb_dl = na_change_type(rb_dl, NA_SFLOAT);
  dl = NA_PTR_TYPE(rb_dl, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 2);
  ldx = NA_SHAPE0(rb_x);
  nrhs = NA_SHAPE1(rb_x);
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (7th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SFLOAT)
    rb_du = na_change_type(rb_du, NA_SFLOAT);
  du = NA_PTR_TYPE(rb_du, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
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

  slagtm_(&trans, &n, &nrhs, &alpha, dl, d, du, x, &ldx, &beta, b, &ldb);

  return rb_b;
}

void
init_lapack_slagtm(VALUE mLapack){
  rb_define_module_function(mLapack, "slagtm", rb_slagtm, -1);
}
