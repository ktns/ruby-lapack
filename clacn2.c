#include "rb_lapack.h"

static VALUE
rb_clacn2(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  complex *x; 
  VALUE rb_est;
  real est; 
  VALUE rb_kase;
  integer kase; 
  VALUE rb_isave;
  integer *isave; 
  VALUE rb_x_out__;
  complex *x_out__;
  VALUE rb_isave_out__;
  integer *isave_out__;
  complex *v;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, est, kase, isave = NumRu::Lapack.clacn2( x, est, kase, isave)\n    or\n  NumRu::Lapack.clacn2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLACN2( N, V, X, EST, KASE, ISAVE )\n\n*  Purpose\n*  =======\n*\n*  CLACN2 estimates the 1-norm of a square, complex matrix A.\n*  Reverse communication is used for evaluating matrix-vector products.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The order of the matrix.  N >= 1.\n*\n*  V      (workspace) COMPLEX array, dimension (N)\n*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)\n*         (W is not returned).\n*\n*  X      (input/output) COMPLEX array, dimension (N)\n*         On an intermediate return, X should be overwritten by\n*               A * X,   if KASE=1,\n*               A' * X,  if KASE=2,\n*         where A' is the conjugate transpose of A, and CLACN2 must be\n*         re-called with all the other parameters unchanged.\n*\n*  EST    (input/output) REAL\n*         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be\n*         unchanged from the previous call to CLACN2.\n*         On exit, EST is an estimate (a lower bound) for norm(A). \n*\n*  KASE   (input/output) INTEGER\n*         On the initial call to CLACN2, KASE should be 0.\n*         On an intermediate return, KASE will be 1 or 2, indicating\n*         whether X should be overwritten by A * X  or A' * X.\n*         On the final return from CLACN2, KASE will again be 0.\n*\n*  ISAVE  (input/output) INTEGER array, dimension (3)\n*         ISAVE is used to save variables between calls to SLACN2\n*\n\n*  Further Details\n*  ======= =======\n*\n*  Contributed by Nick Higham, University of Manchester.\n*  Originally named CONEST, dated March 16, 1988.\n*\n*  Reference: N.J. Higham, \"FORTRAN codes for estimating the one-norm of\n*  a real or complex matrix, with applications to condition estimation\",\n*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.\n*\n*  Last modified:  April, 1999\n*\n*  This is a thread safe version of CLACON, which uses the array ISAVE\n*  in place of a SAVE statement, as follows:\n*\n*     CLACON     CLACN2\n*      JUMP     ISAVE(1)\n*      J        ISAVE(2)\n*      ITER     ISAVE(3)\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_x = argv[0];
  rb_est = argv[1];
  rb_kase = argv[2];
  rb_isave = argv[3];

  est = (real)NUM2DBL(rb_est);
  kase = NUM2INT(rb_kase);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  if (!NA_IsNArray(rb_isave))
    rb_raise(rb_eArgError, "isave (4th argument) must be NArray");
  if (NA_RANK(rb_isave) != 1)
    rb_raise(rb_eArgError, "rank of isave (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_isave) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of isave must be %d", 3);
  if (NA_TYPE(rb_isave) != NA_LINT)
    rb_isave = na_change_type(rb_isave, NA_LINT);
  isave = NA_PTR_TYPE(rb_isave, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_x_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, complex*);
  MEMCPY(x_out__, x, complex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = DIM_LEN(3);
    rb_isave_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isave_out__ = NA_PTR_TYPE(rb_isave_out__, integer*);
  MEMCPY(isave_out__, isave, integer, NA_TOTAL(rb_isave));
  rb_isave = rb_isave_out__;
  isave = isave_out__;
  v = ALLOC_N(complex, (n));

  clacn2_(&n, v, x, &est, &kase, isave);

  free(v);
  rb_est = rb_float_new((double)est);
  rb_kase = INT2NUM(kase);
  return rb_ary_new3(4, rb_x, rb_est, rb_kase, rb_isave);
}

void
init_lapack_clacn2(VALUE mLapack){
  rb_define_module_function(mLapack, "clacn2", rb_clacn2, -1);
}
