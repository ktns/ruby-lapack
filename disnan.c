#include "rb_lapack.h"

extern logical disnan_(doublereal *din);

static VALUE
rb_disnan(int argc, VALUE *argv, VALUE self){
  VALUE rb_din;
  doublereal din; 
  VALUE rb___out__;
  logical __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.disnan( din)\n    or\n  NumRu::Lapack.disnan  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_din = argv[0];

  din = NUM2DBL(rb_din);

  __out__ = disnan_(&din);

  rb___out__ = __out__ ? Qtrue : Qfalse;
  return rb___out__;
}

void
init_lapack_disnan(VALUE mLapack){
  rb_define_module_function(mLapack, "disnan", rb_disnan, -1);
}
