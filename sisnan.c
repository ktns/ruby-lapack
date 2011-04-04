#include "rb_lapack.h"

extern logical sisnan_(real *sin);

static VALUE
rb_sisnan(int argc, VALUE *argv, VALUE self){
  VALUE rb_sin;
  real sin; 
  VALUE rb___out__;
  logical __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.sisnan( sin)\n    or\n  NumRu::Lapack.sisnan  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_sin = argv[0];

  sin = (real)NUM2DBL(rb_sin);

  __out__ = sisnan_(&sin);

  rb___out__ = __out__ ? Qtrue : Qfalse;
  return rb___out__;
}

void
init_lapack_sisnan(VALUE mLapack){
  rb_define_module_function(mLapack, "sisnan", rb_sisnan, -1);
}
