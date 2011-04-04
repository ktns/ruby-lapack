#include "rb_lapack.h"

extern integer ilauplo_(char *uplo);

static VALUE
rb_ilauplo(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilauplo( uplo)\n    or\n  NumRu::Lapack.ilauplo  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_uplo = argv[0];

  uplo = StringValueCStr(rb_uplo)[0];

  __out__ = ilauplo_(&uplo);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_ilauplo(VALUE mLapack){
  rb_define_module_function(mLapack, "ilauplo", rb_ilauplo, -1);
}
