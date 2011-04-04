#include "rb_lapack.h"

extern real slapy3_(real *x, real *y, real *z);

static VALUE
rb_slapy3(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  real x; 
  VALUE rb_y;
  real y; 
  VALUE rb_z;
  real z; 
  VALUE rb___out__;
  real __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slapy3( x, y, z)\n    or\n  NumRu::Lapack.slapy3  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_x = argv[0];
  rb_y = argv[1];
  rb_z = argv[2];

  x = (real)NUM2DBL(rb_x);
  y = (real)NUM2DBL(rb_y);
  z = (real)NUM2DBL(rb_z);

  __out__ = slapy3_(&x, &y, &z);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_slapy3(VALUE mLapack){
  rb_define_module_function(mLapack, "slapy3", rb_slapy3, -1);
}
