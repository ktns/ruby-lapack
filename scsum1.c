#include "rb_lapack.h"

extern real scsum1_(integer *n, complex *cx, integer *incx);

static VALUE
rb_scsum1(int argc, VALUE *argv, VALUE self){
  VALUE rb_cx;
  complex *cx; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb___out__;
  real __out__; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.scsum1( cx, incx)\n    or\n  NumRu::Lapack.scsum1  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_cx = argv[0];
  rb_incx = argv[1];

  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_cx))
    rb_raise(rb_eArgError, "cx (1th argument) must be NArray");
  if (NA_RANK(rb_cx) != 1)
    rb_raise(rb_eArgError, "rank of cx (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_cx);
  if (NA_TYPE(rb_cx) != NA_SCOMPLEX)
    rb_cx = na_change_type(rb_cx, NA_SCOMPLEX);
  cx = NA_PTR_TYPE(rb_cx, complex*);

  __out__ = scsum1_(&n, cx, &incx);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_scsum1(VALUE mLapack){
  rb_define_module_function(mLapack, "scsum1", rb_scsum1, -1);
}
