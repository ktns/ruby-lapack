#include "rb_lapack.h"

extern integer iparmq_(integer *ispec, char *name, char *opts, integer *n, integer *ilo, integer *ihi, integer *lwork);

static VALUE
rb_iparmq(int argc, VALUE *argv, VALUE self){
  VALUE rb_ispec;
  integer ispec; 
  VALUE rb_name;
  char name; 
  VALUE rb_opts;
  char opts; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.iparmq( ispec, name, opts, n, ilo, ihi, lwork)\n    or\n  NumRu::Lapack.iparmq  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_ispec = argv[0];
  rb_name = argv[1];
  rb_opts = argv[2];
  rb_n = argv[3];
  rb_ilo = argv[4];
  rb_ihi = argv[5];
  rb_lwork = argv[6];

  name = StringValueCStr(rb_name)[0];
  ilo = NUM2INT(rb_ilo);
  opts = StringValueCStr(rb_opts)[0];
  n = NUM2INT(rb_n);
  lwork = NUM2INT(rb_lwork);
  ispec = NUM2INT(rb_ispec);
  ihi = NUM2INT(rb_ihi);

  __out__ = iparmq_(&ispec, &name, &opts, &n, &ilo, &ihi, &lwork);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_iparmq(VALUE mLapack){
  rb_define_module_function(mLapack, "iparmq", rb_iparmq, -1);
}
