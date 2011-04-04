#include "rb_lapack.h"

extern doublereal dlapy2_(doublereal *x, doublereal *y);

static VALUE
rb_dlapy2(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  doublereal x; 
  VALUE rb_y;
  doublereal y; 
  VALUE rb___out__;
  doublereal __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlapy2( x, y)\n    or\n  NumRu::Lapack.dlapy2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_x = argv[0];
  rb_y = argv[1];

  x = NUM2DBL(rb_x);
  y = NUM2DBL(rb_y);

  __out__ = dlapy2_(&x, &y);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dlapy2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlapy2", rb_dlapy2, -1);
}
