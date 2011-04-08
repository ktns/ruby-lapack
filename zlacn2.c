#include "rb_lapack.h"

extern VOID zlacn2_(integer *n, doublecomplex *v, doublecomplex *x, doublereal *est, integer *kase, integer *isave);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlacn2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  doublecomplex *x; 
  VALUE rblapack_est;
  doublereal est; 
  VALUE rblapack_kase;
  integer kase; 
  VALUE rblapack_isave;
  integer *isave; 
  VALUE rblapack_x_out__;
  doublecomplex *x_out__;
  VALUE rblapack_isave_out__;
  integer *isave_out__;
  doublecomplex *v;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, est, kase, isave = NumRu::Lapack.zlacn2( x, est, kase, isave, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLACN2( N, V, X, EST, KASE, ISAVE )\n\n*  Purpose\n*  =======\n*\n*  ZLACN2 estimates the 1-norm of a square, complex matrix A.\n*  Reverse communication is used for evaluating matrix-vector products.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The order of the matrix.  N >= 1.\n*\n*  V      (workspace) COMPLEX*16 array, dimension (N)\n*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)\n*         (W is not returned).\n*\n*  X      (input/output) COMPLEX*16 array, dimension (N)\n*         On an intermediate return, X should be overwritten by\n*               A * X,   if KASE=1,\n*               A' * X,  if KASE=2,\n*         where A' is the conjugate transpose of A, and ZLACN2 must be\n*         re-called with all the other parameters unchanged.\n*\n*  EST    (input/output) DOUBLE PRECISION\n*         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be\n*         unchanged from the previous call to ZLACN2.\n*         On exit, EST is an estimate (a lower bound) for norm(A). \n*\n*  KASE   (input/output) INTEGER\n*         On the initial call to ZLACN2, KASE should be 0.\n*         On an intermediate return, KASE will be 1 or 2, indicating\n*         whether X should be overwritten by A * X  or A' * X.\n*         On the final return from ZLACN2, KASE will again be 0.\n*\n*  ISAVE  (input/output) INTEGER array, dimension (3)\n*         ISAVE is used to save variables between calls to ZLACN2\n*\n\n*  Further Details\n*  ======= =======\n*\n*  Contributed by Nick Higham, University of Manchester.\n*  Originally named CONEST, dated March 16, 1988.\n*\n*  Reference: N.J. Higham, \"FORTRAN codes for estimating the one-norm of\n*  a real or complex matrix, with applications to condition estimation\",\n*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.\n*\n*  Last modified:  April, 1999\n*\n*  This is a thread safe version of ZLACON, which uses the array ISAVE\n*  in place of a SAVE statement, as follows:\n*\n*     ZLACON     ZLACN2\n*      JUMP     ISAVE(1)\n*      J        ISAVE(2)\n*      ITER     ISAVE(3)\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, est, kase, isave = NumRu::Lapack.zlacn2( x, est, kase, isave, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_x = argv[0];
  rblapack_est = argv[1];
  rblapack_kase = argv[2];
  rblapack_isave = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_DCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  est = NUM2DBL(rblapack_est);
  if (!NA_IsNArray(rblapack_isave))
    rb_raise(rb_eArgError, "isave (4th argument) must be NArray");
  if (NA_RANK(rblapack_isave) != 1)
    rb_raise(rb_eArgError, "rank of isave (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_isave) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of isave must be %d", 3);
  if (NA_TYPE(rblapack_isave) != NA_LINT)
    rblapack_isave = na_change_type(rblapack_isave, NA_LINT);
  isave = NA_PTR_TYPE(rblapack_isave, integer*);
  kase = NUM2INT(rblapack_kase);
  {
    int shape[1];
    shape[0] = n;
    rblapack_x_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 3;
    rblapack_isave_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isave_out__ = NA_PTR_TYPE(rblapack_isave_out__, integer*);
  MEMCPY(isave_out__, isave, integer, NA_TOTAL(rblapack_isave));
  rblapack_isave = rblapack_isave_out__;
  isave = isave_out__;
  v = ALLOC_N(doublecomplex, (n));

  zlacn2_(&n, v, x, &est, &kase, isave);

  free(v);
  rblapack_est = rb_float_new((double)est);
  rblapack_kase = INT2NUM(kase);
  return rb_ary_new3(4, rblapack_x, rblapack_est, rblapack_kase, rblapack_isave);
}

void
init_lapack_zlacn2(VALUE mLapack){
  rb_define_module_function(mLapack, "zlacn2", rblapack_zlacn2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
