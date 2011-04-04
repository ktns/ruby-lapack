#include "rb_lapack.h"

extern VOID slargv_(integer *n, real *x, integer *incx, real *y, integer *incy, real *c, integer *incc);

static VALUE
rb_slargv(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_x;
  real *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_y;
  real *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_incc;
  integer incc; 
  VALUE rb_c;
  real *c; 
  VALUE rb_x_out__;
  real *x_out__;
  VALUE rb_y_out__;
  real *y_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c, x, y = NumRu::Lapack.slargv( n, x, incx, y, incy, incc)\n    or\n  NumRu::Lapack.slargv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_n = argv[0];
  rb_x = argv[1];
  rb_incx = argv[2];
  rb_y = argv[3];
  rb_incy = argv[4];
  rb_incc = argv[5];

  incy = NUM2INT(rb_incy);
  incc = NUM2INT(rb_incc);
  n = NUM2INT(rb_n);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (4th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (1+(n-1)*incy))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incy);
  if (NA_TYPE(rb_y) != NA_SFLOAT)
    rb_y = na_change_type(rb_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rb_y, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incc;
    rb_c = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, real*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rb_x_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incy;
    rb_y_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, real*);
  MEMCPY(y_out__, y, real, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  slargv_(&n, x, &incx, y, &incy, c, &incc);

  return rb_ary_new3(3, rb_c, rb_x, rb_y);
}

void
init_lapack_slargv(VALUE mLapack){
  rb_define_module_function(mLapack, "slargv", rb_slargv, -1);
}
