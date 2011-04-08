#include "rb_lapack.h"

extern VOID cladiv_(complex *__out__, complex *x, complex *y);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cladiv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  complex x; 
  VALUE rblapack_y;
  complex y; 
  VALUE rblapack___out__;
  complex __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.cladiv( x, y, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      COMPLEX FUNCTION CLADIV( X, Y )\n\n*  Purpose\n*  =======\n*\n*  CLADIV := X / Y, where X and Y are complex.  The computation of X / Y\n*  will not overflow on an intermediary step unless the results\n*  overflows.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) COMPLEX\n*  Y       (input) COMPLEX\n*          The complex scalars X and Y.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      REAL               ZI, ZR\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SLADIV\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          AIMAG, CMPLX, REAL\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.cladiv( x, y, [:usage => usage, :help => help])\n");
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

  x.r = (real)NUM2DBL(rb_funcall(rblapack_x, rb_intern("real"), 0));
  x.i = (real)NUM2DBL(rb_funcall(rblapack_x, rb_intern("imag"), 0));
  y.r = (real)NUM2DBL(rb_funcall(rblapack_y, rb_intern("real"), 0));
  y.i = (real)NUM2DBL(rb_funcall(rblapack_y, rb_intern("imag"), 0));

  cladiv_(&__out__, &x, &y);

  rblapack___out__ = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(__out__.r)), rb_float_new((double)(__out__.i)));
  return rblapack___out__;
}

void
init_lapack_cladiv(VALUE mLapack){
  rb_define_module_function(mLapack, "cladiv", rblapack_cladiv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
