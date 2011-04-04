#include "rb_lapack.h"

extern logical lsamen_(integer *n, char *ca, char *cb);

static VALUE
rb_lsamen(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_ca;
  char *ca; 
  VALUE rb_cb;
  char *cb; 
  VALUE rb___out__;
  logical __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.lsamen( n, ca, cb)\n    or\n  NumRu::Lapack.lsamen  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_n = argv[0];
  rb_ca = argv[1];
  rb_cb = argv[2];

  n = NUM2INT(rb_n);
  ca = StringValueCStr(rb_ca);
  cb = StringValueCStr(rb_cb);

  __out__ = lsamen_(&n, ca, cb);

  rb___out__ = __out__ ? Qtrue : Qfalse;
  return rb___out__;
}

void
init_lapack_lsamen(VALUE mLapack){
  rb_define_module_function(mLapack, "lsamen", rb_lsamen, -1);
}
