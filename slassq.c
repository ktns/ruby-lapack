#include "rb_lapack.h"

extern VOID slassq_(integer *n, real *x, integer *incx, real *scale, real *sumsq);

static VALUE
rb_slassq(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  real *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_scale;
  real scale; 
  VALUE rb_sumsq;
  real sumsq; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, sumsq = NumRu::Lapack.slassq( x, incx, scale, sumsq)\n    or\n  NumRu::Lapack.slassq  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_x = argv[0];
  rb_incx = argv[1];
  rb_scale = argv[2];
  rb_sumsq = argv[3];

  scale = (real)NUM2DBL(rb_scale);
  sumsq = (real)NUM2DBL(rb_sumsq);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  incx = NUM2INT(rb_incx);

  slassq_(&n, x, &incx, &scale, &sumsq);

  rb_scale = rb_float_new((double)scale);
  rb_sumsq = rb_float_new((double)sumsq);
  return rb_ary_new3(2, rb_scale, rb_sumsq);
}

void
init_lapack_slassq(VALUE mLapack){
  rb_define_module_function(mLapack, "slassq", rb_slassq, -1);
}
