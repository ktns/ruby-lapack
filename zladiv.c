#include "rb_lapack.h"

extern VOID zladiv_(doublecomplex *__out__, doublecomplex *x, doublecomplex *y);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zladiv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  doublecomplex x; 
  VALUE rblapack_y;
  doublecomplex y; 
  VALUE rblapack___out__;
  doublecomplex __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zladiv( x, y, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      COMPLEX*16     FUNCTION ZLADIV( X, Y )\n\n*  Purpose\n*  =======\n*\n*  ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y\n*  will not overflow on an intermediary step unless the results\n*  overflows.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) COMPLEX*16\n*  Y       (input) COMPLEX*16\n*          The complex scalars X and Y.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      DOUBLE PRECISION   ZI, ZR\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           DLADIV\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          DBLE, DCMPLX, DIMAG\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zladiv( x, y, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_x = argv[0];
  rblapack_y = argv[1];
  if (rb_options != Qnil) {
  }

  x.r = NUM2DBL(rb_funcall(rblapack_x, rb_intern("real"), 0));
  x.i = NUM2DBL(rb_funcall(rblapack_x, rb_intern("imag"), 0));
  y.r = NUM2DBL(rb_funcall(rblapack_y, rb_intern("real"), 0));
  y.i = NUM2DBL(rb_funcall(rblapack_y, rb_intern("imag"), 0));

  zladiv_(&__out__, &x, &y);

  rblapack___out__ = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(__out__.r)), rb_float_new((double)(__out__.i)));
  return rblapack___out__;
}

void
init_lapack_zladiv(VALUE mLapack){
  rb_define_module_function(mLapack, "zladiv", rblapack_zladiv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
