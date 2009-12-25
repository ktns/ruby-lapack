#include "rb_lapack.h"

static VALUE
rb_cla_wwaddw(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  complex *x; 
  VALUE rb_y;
  complex *y; 
  VALUE rb_w;
  complex *w; 
  VALUE rb_x_out__;
  complex *x_out__;
  VALUE rb_y_out__;
  complex *y_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, y = NumRu::Lapack.cla_wwaddw( x, y, w)\n    or\n  NumRu::Lapack.cla_wwaddw  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLA_WWADDW( N, X, Y, W )\n\n*     Purpose\n*     =======\n*\n*     CLA_WWADDW adds a vector W into a doubled-single vector (X, Y).\n*\n*     This works for all extant IBM's hex and binary floating point\n*     arithmetics, but not for decimal.\n*\n\n*     Arguments\n*     =========\n*\n*     N      (input) INTEGER\n*            The length of vectors X, Y, and W.\n*\n*     X, Y   (input/output) COMPLEX array, length N\n*            The doubled-single accumulation vector.\n*\n*     W      (input) COMPLEX array, length N\n*            The vector to be added.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      COMPLEX            S\n      INTEGER            I\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_x = argv[0];
  rb_y = argv[1];
  rb_w = argv[2];

  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (2th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of y must be the same as shape 0 of x");
  if (NA_TYPE(rb_y) != NA_SCOMPLEX)
    rb_y = na_change_type(rb_y, NA_SCOMPLEX);
  y = NA_PTR_TYPE(rb_y, complex*);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (3th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_w) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of w must be the same as shape 0 of x");
  if (NA_TYPE(rb_w) != NA_SCOMPLEX)
    rb_w = na_change_type(rb_w, NA_SCOMPLEX);
  w = NA_PTR_TYPE(rb_w, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_x_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, complex*);
  MEMCPY(x_out__, x, complex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_y_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, complex*);
  MEMCPY(y_out__, y, complex, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  cla_wwaddw_(&n, x, y, w);

  return rb_ary_new3(2, rb_x, rb_y);
}

void
init_lapack_cla_wwaddw(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_wwaddw", rb_cla_wwaddw, -1);
}
