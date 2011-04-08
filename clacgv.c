#include "rb_lapack.h"

extern VOID clacgv_(integer *n, complex *x, integer *incx);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clacgv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_x;
  complex *x; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_x_out__;
  complex *x_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x = NumRu::Lapack.clacgv( n, x, incx, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLACGV( N, X, INCX )\n\n*  Purpose\n*  =======\n*\n*  CLACGV conjugates a complex vector of length N.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The length of the vector X.  N >= 0.\n*\n*  X       (input/output) COMPLEX array, dimension\n*                         (1+(N-1)*abs(INCX))\n*          On entry, the vector of length N to be conjugated.\n*          On exit, X is overwritten with conjg(X).\n*\n*  INCX    (input) INTEGER\n*          The spacing between successive elements of X.\n*\n\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IOFF\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          CONJG\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x = NumRu::Lapack.clacgv( n, x, incx, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_n = argv[0];
  rblapack_x = argv[1];
  rblapack_incx = argv[2];
  if (rb_options != Qnil) {
  }

  incx = NUM2INT(rblapack_incx);
  n = NUM2INT(rblapack_n);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (1+(n-1)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*abs(incx));
  if (NA_TYPE(rblapack_x) != NA_SCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, complex*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*abs(incx);
    rblapack_x_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, complex*);
  MEMCPY(x_out__, x, complex, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;

  clacgv_(&n, x, &incx);

  return rblapack_x;
}

void
init_lapack_clacgv(VALUE mLapack){
  rb_define_module_function(mLapack, "clacgv", rblapack_clacgv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
