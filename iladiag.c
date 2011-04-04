#include "rb_lapack.h"

extern integer iladiag_(char *diag);

static VALUE
rb_iladiag(int argc, VALUE *argv, VALUE self){
  VALUE rb_diag;
  char diag; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.iladiag( diag)\n    or\n  NumRu::Lapack.iladiag  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_diag = argv[0];

  diag = StringValueCStr(rb_diag)[0];

  __out__ = iladiag_(&diag);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_iladiag(VALUE mLapack){
  rb_define_module_function(mLapack, "iladiag", rb_iladiag, -1);
}
