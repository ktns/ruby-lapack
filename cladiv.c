#include "rb_lapack.h"

extern VOID cladiv_(complex *__out__, complex *x, complex *y);

static VALUE
rb_cladiv(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  complex x; 
  VALUE rb_y;
  complex y; 
  VALUE rb___out__;
  complex __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.cladiv( x, y)\n    or\n  NumRu::Lapack.cladiv  # print help\n\n\nFORTRAN MANUAL\n      COMPLEX FUNCTION CLADIV( X, Y )\n\n*  Purpose\n*  =======\n*\n*  CLADIV := X / Y, where X and Y are complex.  The computation of X / Y\n*  will not overflow on an intermediary step unless the results\n*  overflows.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) COMPLEX\n*  Y       (input) COMPLEX\n*          The complex scalars X and Y.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      REAL               ZI, ZR\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SLADIV\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          AIMAG, CMPLX, REAL\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_x = argv[0];
  rb_y = argv[1];

  x.r = (real)NUM2DBL(rb_funcall(rb_x, rb_intern("real"), 0));
  x.i = (real)NUM2DBL(rb_funcall(rb_x, rb_intern("imag"), 0));
  y.r = (real)NUM2DBL(rb_funcall(rb_y, rb_intern("real"), 0));
  y.i = (real)NUM2DBL(rb_funcall(rb_y, rb_intern("imag"), 0));

  cladiv_(&__out__, &x, &y);

  rb___out__ = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(__out__.r)), rb_float_new((double)(__out__.i)));
  return rb___out__;
}

void
init_lapack_cladiv(VALUE mLapack){
  rb_define_module_function(mLapack, "cladiv", rb_cladiv, -1);
}
