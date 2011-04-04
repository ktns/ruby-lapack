#include "rb_lapack.h"

extern VOID zlacon_(integer *n, doublecomplex *v, doublecomplex *x, doublereal *est, integer *kase);

static VALUE
rb_zlacon(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_est;
  doublereal est; 
  VALUE rb_kase;
  integer kase; 
  VALUE rb_x_out__;
  doublecomplex *x_out__;
  doublecomplex *v;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, est, kase = NumRu::Lapack.zlacon( x, est, kase)\n    or\n  NumRu::Lapack.zlacon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLACON( N, V, X, EST, KASE )\n\n*  Purpose\n*  =======\n*\n*  ZLACON estimates the 1-norm of a square, complex matrix A.\n*  Reverse communication is used for evaluating matrix-vector products.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The order of the matrix.  N >= 1.\n*\n*  V      (workspace) COMPLEX*16 array, dimension (N)\n*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)\n*         (W is not returned).\n*\n*  X      (input/output) COMPLEX*16 array, dimension (N)\n*         On an intermediate return, X should be overwritten by\n*               A * X,   if KASE=1,\n*               A' * X,  if KASE=2,\n*         where A' is the conjugate transpose of A, and ZLACON must be\n*         re-called with all the other parameters unchanged.\n*\n*  EST    (input/output) DOUBLE PRECISION\n*         On entry with KASE = 1 or 2 and JUMP = 3, EST should be\n*         unchanged from the previous call to ZLACON.\n*         On exit, EST is an estimate (a lower bound) for norm(A). \n*\n*  KASE   (input/output) INTEGER\n*         On the initial call to ZLACON, KASE should be 0.\n*         On an intermediate return, KASE will be 1 or 2, indicating\n*         whether X should be overwritten by A * X  or A' * X.\n*         On the final return from ZLACON, KASE will again be 0.\n*\n\n*  Further Details\n*  ======= =======\n*\n*  Contributed by Nick Higham, University of Manchester.\n*  Originally named CONEST, dated March 16, 1988.\n*\n*  Reference: N.J. Higham, \"FORTRAN codes for estimating the one-norm of\n*  a real or complex matrix, with applications to condition estimation\",\n*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.\n*\n*  Last modified:  April, 1999\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_x = argv[0];
  rb_est = argv[1];
  rb_kase = argv[2];

  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  est = NUM2DBL(rb_est);
  kase = NUM2INT(rb_kase);
  {
    int shape[1];
    shape[0] = n;
    rb_x_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  v = ALLOC_N(doublecomplex, (n));

  zlacon_(&n, v, x, &est, &kase);

  free(v);
  rb_est = rb_float_new((double)est);
  rb_kase = INT2NUM(kase);
  return rb_ary_new3(3, rb_x, rb_est, rb_kase);
}

void
init_lapack_zlacon(VALUE mLapack){
  rb_define_module_function(mLapack, "zlacon", rb_zlacon, -1);
}
