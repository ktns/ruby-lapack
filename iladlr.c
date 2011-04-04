#include "rb_lapack.h"

extern integer iladlr_(integer *m, integer *n, doublereal *a, integer *lda);

static VALUE
rb_iladlr(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb___out__;
  integer __out__; 

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.iladlr( m, a)\n    or\n  NumRu::Lapack.iladlr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_m = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  m = NUM2INT(rb_m);

  __out__ = iladlr_(&m, &n, a, &lda);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_iladlr(VALUE mLapack){
  rb_define_module_function(mLapack, "iladlr", rb_iladlr, -1);
}
