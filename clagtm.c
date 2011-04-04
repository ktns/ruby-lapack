#include "rb_lapack.h"

extern VOID clagtm_(char *trans, integer *n, integer *nrhs, real *alpha, complex *dl, complex *d, complex *du, complex *x, integer *ldx, real *beta, complex *b, integer *ldb);

static VALUE
rb_clagtm(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_alpha;
  real alpha; 
  VALUE rb_dl;
  complex *dl; 
  VALUE rb_d;
  complex *d; 
  VALUE rb_du;
  complex *du; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_beta;
  real beta; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_b_out__;
  complex *b_out__;

  integer n;
  integer ldx;
  integer nrhs;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  b = NumRu::Lapack.clagtm( trans, alpha, dl, d, du, x, beta, b)\n    or\n  NumRu::Lapack.clagtm  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, B, LDB )\n\n*  Purpose\n*  =======\n*\n*  CLAGTM performs a matrix-vector product of the form\n*\n*     B := alpha * A * X + beta * B\n*\n*  where A is a tridiagonal matrix of order N, B and X are N by NRHS\n*  matrices, and alpha and beta are real scalars, each of which may be\n*  0., 1., or -1.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the operation applied to A.\n*          = 'N':  No transpose, B := alpha * A * X + beta * B\n*          = 'T':  Transpose,    B := alpha * A**T * X + beta * B\n*          = 'C':  Conjugate transpose, B := alpha * A**H * X + beta * B\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices X and B.\n*\n*  ALPHA   (input) REAL\n*          The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise,\n*          it is assumed to be 0.\n*\n*  DL      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) sub-diagonal elements of T.\n*\n*  D       (input) COMPLEX array, dimension (N)\n*          The diagonal elements of T.\n*\n*  DU      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) super-diagonal elements of T.\n*\n*  X       (input) COMPLEX array, dimension (LDX,NRHS)\n*          The N by NRHS matrix X.\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(N,1).\n*\n*  BETA    (input) REAL\n*          The scalar beta.  BETA must be 0., 1., or -1.; otherwise,\n*          it is assumed to be 1.\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)\n*          On entry, the N by NRHS matrix B.\n*          On exit, B is overwritten by the matrix expression\n*          B := alpha * A * X + beta * B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(N,1).\n*\n\n*  =====================================================================\n*\n\n");
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

  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of x must be the same as shape 1 of b");
  ldx = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SCOMPLEX)
    rb_d = na_change_type(rb_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rb_d, complex*);
  beta = (real)NUM2DBL(rb_beta);
  alpha = (real)NUM2DBL(rb_alpha);
  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (3th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SCOMPLEX)
    rb_dl = na_change_type(rb_dl, NA_SCOMPLEX);
  dl = NA_PTR_TYPE(rb_dl, complex*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (5th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SCOMPLEX)
    rb_du = na_change_type(rb_du, NA_SCOMPLEX);
  du = NA_PTR_TYPE(rb_du, complex*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  clagtm_(&trans, &n, &nrhs, &alpha, dl, d, du, x, &ldx, &beta, b, &ldb);

  return rb_b;
}

void
init_lapack_clagtm(VALUE mLapack){
  rb_define_module_function(mLapack, "clagtm", rb_clagtm, -1);
}
