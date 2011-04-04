#include "rb_lapack.h"

extern real slapy2_(real *x, real *y);

static VALUE
rb_slapy2(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  real x; 
  VALUE rb_y;
  real y; 
  VALUE rb___out__;
  real __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slapy2( x, y)\n    or\n  NumRu::Lapack.slapy2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_x = argv[0];
  rb_y = argv[1];

  x = (real)NUM2DBL(rb_x);
  y = (real)NUM2DBL(rb_y);

  __out__ = slapy2_(&x, &y);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_slapy2(VALUE mLapack){
  rb_define_module_function(mLapack, "slapy2", rb_slapy2, -1);
}
