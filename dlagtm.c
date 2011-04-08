#include "rb_lapack.h"

extern VOID dlagtm_(char *trans, integer *n, integer *nrhs, doublereal *alpha, doublereal *dl, doublereal *d, doublereal *du, doublereal *x, integer *ldx, doublereal *beta, doublereal *b, integer *ldb);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlagtm(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_alpha;
  doublereal alpha; 
  VALUE rblapack_dl;
  doublereal *dl; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_du;
  doublereal *du; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_beta;
  doublereal beta; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_b_out__;
  doublereal *b_out__;

  integer n;
  integer ldx;
  integer nrhs;
  integer ldb;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  b = NumRu::Lapack.dlagtm( trans, alpha, dl, d, du, x, beta, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, B, LDB )\n\n*  Purpose\n*  =======\n*\n*  DLAGTM performs a matrix-vector product of the form\n*\n*     B := alpha * A * X + beta * B\n*\n*  where A is a tridiagonal matrix of order N, B and X are N by NRHS\n*  matrices, and alpha and beta are real scalars, each of which may be\n*  0., 1., or -1.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the operation applied to A.\n*          = 'N':  No transpose, B := alpha * A * X + beta * B\n*          = 'T':  Transpose,    B := alpha * A'* X + beta * B\n*          = 'C':  Conjugate transpose = Transpose\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices X and B.\n*\n*  ALPHA   (input) DOUBLE PRECISION\n*          The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise,\n*          it is assumed to be 0.\n*\n*  DL      (input) DOUBLE PRECISION array, dimension (N-1)\n*          The (n-1) sub-diagonal elements of T.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The diagonal elements of T.\n*\n*  DU      (input) DOUBLE PRECISION array, dimension (N-1)\n*          The (n-1) super-diagonal elements of T.\n*\n*  X       (input) DOUBLE PRECISION array, dimension (LDX,NRHS)\n*          The N by NRHS matrix X.\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(N,1).\n*\n*  BETA    (input) DOUBLE PRECISION\n*          The scalar beta.  BETA must be 0., 1., or -1.; otherwise,\n*          it is assumed to be 1.\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          On entry, the N by NRHS matrix B.\n*          On exit, B is overwritten by the matrix expression\n*          B := alpha * A * X + beta * B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(N,1).\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  b = NumRu::Lapack.dlagtm( trans, alpha, dl, d, du, x, beta, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_trans = argv[0];
  rblapack_alpha = argv[1];
  rblapack_dl = argv[2];
  rblapack_d = argv[3];
  rblapack_du = argv[4];
  rblapack_x = argv[5];
  rblapack_beta = argv[6];
  rblapack_b = argv[7];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 2)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_x) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of x must be the same as shape 1 of b");
  ldx = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_DFLOAT)
    rblapack_x = na_change_type(rblapack_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rblapack_x, doublereal*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  beta = NUM2DBL(rblapack_beta);
  alpha = NUM2DBL(rblapack_alpha);
  trans = StringValueCStr(rblapack_trans)[0];
  if (!NA_IsNArray(rblapack_dl))
    rb_raise(rb_eArgError, "dl (3th argument) must be NArray");
  if (NA_RANK(rblapack_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rblapack_dl) != NA_DFLOAT)
    rblapack_dl = na_change_type(rblapack_dl, NA_DFLOAT);
  dl = NA_PTR_TYPE(rblapack_dl, doublereal*);
  if (!NA_IsNArray(rblapack_du))
    rb_raise(rb_eArgError, "du (5th argument) must be NArray");
  if (NA_RANK(rblapack_du) != 1)
    rb_raise(rb_eArgError, "rank of du (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rblapack_du) != NA_DFLOAT)
    rblapack_du = na_change_type(rblapack_du, NA_DFLOAT);
  du = NA_PTR_TYPE(rblapack_du, doublereal*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  dlagtm_(&trans, &n, &nrhs, &alpha, dl, d, du, x, &ldx, &beta, b, &ldb);

  return rblapack_b;
}

void
init_lapack_dlagtm(VALUE mLapack){
  rb_define_module_function(mLapack, "dlagtm", rblapack_dlagtm, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
