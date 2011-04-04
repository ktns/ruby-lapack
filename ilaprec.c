#include "rb_lapack.h"

extern integer ilaprec_(char *prec);

static VALUE
rb_ilaprec(int argc, VALUE *argv, VALUE self){
  VALUE rb_prec;
  char prec; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilaprec( prec)\n    or\n  NumRu::Lapack.ilaprec  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_prec = argv[0];

  prec = StringValueCStr(rb_prec)[0];

  __out__ = ilaprec_(&prec);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_ilaprec(VALUE mLapack){
  rb_define_module_function(mLapack, "ilaprec", rb_ilaprec, -1);
}
