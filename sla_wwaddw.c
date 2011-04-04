#include "rb_lapack.h"

extern VOID sla_wwaddw_(integer *n, real *x, real *y, real *w);

static VALUE
rb_sla_wwaddw(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  real *x; 
  VALUE rb_y;
  real *y; 
  VALUE rb_w;
  real *w; 
  VALUE rb_x_out__;
  real *x_out__;
  VALUE rb_y_out__;
  real *y_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, y = NumRu::Lapack.sla_wwaddw( x, y, w)\n    or\n  NumRu::Lapack.sla_wwaddw  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_x = argv[0];
  rb_y = argv[1];
  rb_w = argv[2];

  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (3th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 0 of w");
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (2th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of y must be the same as shape 0 of w");
  if (NA_TYPE(rb_y) != NA_SFLOAT)
    rb_y = na_change_type(rb_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rb_y, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_x_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_y_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, real*);
  MEMCPY(y_out__, y, real, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  sla_wwaddw_(&n, x, y, w);

  return rb_ary_new3(2, rb_x, rb_y);
}

void
init_lapack_sla_wwaddw(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_wwaddw", rb_sla_wwaddw, -1);
}
