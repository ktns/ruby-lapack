#include "rb_lapack.h"

extern integer ieeeck_(integer *ispec, real *zero, real *one);

static VALUE
rb_ieeeck(int argc, VALUE *argv, VALUE self){
  VALUE rb_ispec;
  integer ispec; 
  VALUE rb_zero;
  real zero; 
  VALUE rb_one;
  real one; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ieeeck( ispec, zero, one)\n    or\n  NumRu::Lapack.ieeeck  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_ispec = argv[0];
  rb_zero = argv[1];
  rb_one = argv[2];

  one = (real)NUM2DBL(rb_one);
  ispec = NUM2INT(rb_ispec);
  zero = (real)NUM2DBL(rb_zero);

  __out__ = ieeeck_(&ispec, &zero, &one);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_ieeeck(VALUE mLapack){
  rb_define_module_function(mLapack, "ieeeck", rb_ieeeck, -1);
}
