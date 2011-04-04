#include "rb_lapack.h"

extern VOID dlapll_(integer *n, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *ssmin);

static VALUE
rb_dlapll(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_y;
  doublereal *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_ssmin;
  doublereal ssmin; 
  VALUE rb_x_out__;
  doublereal *x_out__;
  VALUE rb_y_out__;
  doublereal *y_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ssmin, x, y = NumRu::Lapack.dlapll( n, x, incx, y, incy)\n    or\n  NumRu::Lapack.dlapll  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_n = argv[0];
  rb_x = argv[1];
  rb_incx = argv[2];
  rb_y = argv[3];
  rb_incy = argv[4];

  incy = NUM2INT(rb_incy);
  incx = NUM2INT(rb_incx);
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (4th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (1+(n-1)*incy))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incy);
  if (NA_TYPE(rb_y) != NA_DFLOAT)
    rb_y = na_change_type(rb_y, NA_DFLOAT);
  y = NA_PTR_TYPE(rb_y, doublereal*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rb_x) != NA_DFLOAT)
    rb_x = na_change_type(rb_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rb_x_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublereal*);
  MEMCPY(x_out__, x, doublereal, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incy;
    rb_y_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, doublereal*);
  MEMCPY(y_out__, y, doublereal, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  dlapll_(&n, x, &incx, y, &incy, &ssmin);

  rb_ssmin = rb_float_new((double)ssmin);
  return rb_ary_new3(3, rb_ssmin, rb_x, rb_y);
}

void
init_lapack_dlapll(VALUE mLapack){
  rb_define_module_function(mLapack, "dlapll", rb_dlapll, -1);
}
