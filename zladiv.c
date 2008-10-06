#include "rb_lapack.h"

static VALUE
rb_zladiv(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  doublecomplex x; 
  VALUE rb_y;
  doublecomplex y; 
  VALUE rb___out__;
  doublecomplex __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zladiv( x, y)\n    or\n  NumRu::Lapack.zladiv  # print help\n\n\nFORTRAN MANUAL\n      COMPLEX*16     FUNCTION ZLADIV( X, Y )\n\n*  Purpose\n*  =======\n*\n*  ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y\n*  will not overflow on an intermediary step unless the results\n*  overflows.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) COMPLEX*16\n*  Y       (input) COMPLEX*16\n*          The complex scalars X and Y.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      DOUBLE PRECISION   ZI, ZR\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           DLADIV\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          DBLE, DCMPLX, DIMAG\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_x = argv[0];
  rb_y = argv[1];

  x.r = NUM2DBL(rb_funcall(rb_x, rb_intern("real"), 0));
  x.i = NUM2DBL(rb_funcall(rb_x, rb_intern("imag"), 0));
  y.r = NUM2DBL(rb_funcall(rb_y, rb_intern("real"), 0));
  y.i = NUM2DBL(rb_funcall(rb_y, rb_intern("imag"), 0));

  zladiv_(&__out__, &x, &y);

  rb___out__ = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(__out__.r)), rb_float_new((double)(__out__.i)));
  return rb___out__;
}

void
init_lapack_zladiv(VALUE mLapack){
  rb_define_module_function(mLapack, "zladiv", rb_zladiv, -1);
}
