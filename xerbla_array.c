#include "rb_lapack.h"

extern VOID xerbla_array_(char *srname_array, integer *srname_len, integer *info);

static VALUE
rb_xerbla_array(int argc, VALUE *argv, VALUE self){
  VALUE rb_srname_array;
  char *srname_array; 
  VALUE rb_info;
  integer info; 

  integer srname_len;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n   = NumRu::Lapack.xerbla_array( srname_array, info)\n    or\n  NumRu::Lapack.xerbla_array  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_srname_array = argv[0];
  rb_info = argv[1];

  info = NUM2INT(rb_info);
  srname_array = StringValueCStr(rb_srname_array);

  xerbla_array_(srname_array, &srname_len, &info);

  return Qnil;
}

void
init_lapack_xerbla_array(VALUE mLapack){
  rb_define_module_function(mLapack, "xerbla_array", rb_xerbla_array, -1);
}
