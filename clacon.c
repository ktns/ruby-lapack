#include "rb_lapack.h"

extern VOID clacon_(integer *n, complex *v, complex *x, real *est, integer *kase);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clacon(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  complex *x; 
  VALUE rblapack_est;
  real est; 
  VALUE rblapack_kase;
  integer kase; 
  VALUE rblapack_x_out__;
  complex *x_out__;
  complex *v;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, est, kase = NumRu::Lapack.clacon( x, est, kase, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLACON( N, V, X, EST, KASE )\n\n*  Purpose\n*  =======\n*\n*  CLACON estimates the 1-norm of a square, complex matrix A.\n*  Reverse communication is used for evaluating matrix-vector products.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The order of the matrix.  N >= 1.\n*\n*  V      (workspace) COMPLEX array, dimension (N)\n*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)\n*         (W is not returned).\n*\n*  X      (input/output) COMPLEX array, dimension (N)\n*         On an intermediate return, X should be overwritten by\n*               A * X,   if KASE=1,\n*               A' * X,  if KASE=2,\n*         where A' is the conjugate transpose of A, and CLACON must be\n*         re-called with all the other parameters unchanged.\n*\n*  EST    (input/output) REAL\n*         On entry with KASE = 1 or 2 and JUMP = 3, EST should be\n*         unchanged from the previous call to CLACON.\n*         On exit, EST is an estimate (a lower bound) for norm(A). \n*\n*  KASE   (input/output) INTEGER\n*         On the initial call to CLACON, KASE should be 0.\n*         On an intermediate return, KASE will be 1 or 2, indicating\n*         whether X should be overwritten by A * X  or A' * X.\n*         On the final return from CLACON, KASE will again be 0.\n*\n\n*  Further Details\n*  ======= =======\n*\n*  Contributed by Nick Higham, University of Manchester.\n*  Originally named CONEST, dated March 16, 1988.\n*\n*  Reference: N.J. Higham, \"FORTRAN codes for estimating the one-norm of\n*  a real or complex matrix, with applications to condition estimation\",\n*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.\n*\n*  Last modified:  April, 1999\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, est, kase = NumRu::Lapack.clacon( x, est, kase, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_x = argv[0];
  rblapack_est = argv[1];
  rblapack_kase = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_SCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, complex*);
  est = (real)NUM2DBL(rblapack_est);
  kase = NUM2INT(rblapack_kase);
  {
    int shape[1];
    shape[0] = n;
    rblapack_x_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, complex*);
  MEMCPY(x_out__, x, complex, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  v = ALLOC_N(complex, (n));

  clacon_(&n, v, x, &est, &kase);

  free(v);
  rblapack_est = rb_float_new((double)est);
  rblapack_kase = INT2NUM(kase);
  return rb_ary_new3(3, rblapack_x, rblapack_est, rblapack_kase);
}

void
init_lapack_clacon(VALUE mLapack){
  rb_define_module_function(mLapack, "clacon", rblapack_clacon, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
