#include "rb_lapack.h"

extern integer ilaenv_(integer *ispec, char *name, char *opts, integer *n1, integer *n2, integer *n3, integer *n4);

static VALUE
rb_ilaenv(int argc, VALUE *argv, VALUE self){
  VALUE rb_ispec;
  integer ispec; 
  VALUE rb_name;
  char *name; 
  VALUE rb_opts;
  char *opts; 
  VALUE rb_n1;
  integer n1; 
  VALUE rb_n2;
  integer n2; 
  VALUE rb_n3;
  integer n3; 
  VALUE rb_n4;
  integer n4; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilaenv( ispec, name, opts, n1, n2, n3, n4)\n    or\n  NumRu::Lapack.ilaenv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_ispec = argv[0];
  rb_name = argv[1];
  rb_opts = argv[2];
  rb_n1 = argv[3];
  rb_n2 = argv[4];
  rb_n3 = argv[5];
  rb_n4 = argv[6];

  name = StringValueCStr(rb_name);
  opts = StringValueCStr(rb_opts);
  n1 = NUM2INT(rb_n1);
  n2 = NUM2INT(rb_n2);
  n3 = NUM2INT(rb_n3);
  n4 = NUM2INT(rb_n4);
  ispec = NUM2INT(rb_ispec);

  __out__ = ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_ilaenv(VALUE mLapack){
  rb_define_module_function(mLapack, "ilaenv", rb_ilaenv, -1);
}
