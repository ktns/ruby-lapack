#include "rb_lapack.h"

extern VOID clar2v_(integer *n, complex *x, complex *y, complex *z, integer *incx, real *c, complex *s, integer *incc);

static VALUE
rb_clar2v(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_y;
  complex *y; 
  VALUE rb_z;
  complex *z; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_c;
  real *c; 
  VALUE rb_s;
  complex *s; 
  VALUE rb_incc;
  integer incc; 
  VALUE rb_x_out__;
  complex *x_out__;
  VALUE rb_y_out__;
  complex *y_out__;
  VALUE rb_z_out__;
  complex *z_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, y, z = NumRu::Lapack.clar2v( n, x, y, z, incx, c, s, incc)\n    or\n  NumRu::Lapack.clar2v  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_n = argv[0];
  rb_x = argv[1];
  rb_y = argv[2];
  rb_z = argv[3];
  rb_incx = argv[4];
  rb_c = argv[5];
  rb_s = argv[6];
  rb_incc = argv[7];

  n = NUM2INT(rb_n);
  incc = NUM2INT(rb_incc);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (4th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rb_z) != NA_SCOMPLEX)
    rb_z = na_change_type(rb_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rb_z, complex*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (3th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rb_y) != NA_SCOMPLEX)
    rb_y = na_change_type(rb_y, NA_SCOMPLEX);
  y = NA_PTR_TYPE(rb_y, complex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != (1+(n-1)*incc))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", 1+(n-1)*incc);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (7th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != (1+(n-1)*incc))
    rb_raise(rb_eRuntimeError, "shape 0 of s must be %d", 1+(n-1)*incc);
  if (NA_TYPE(rb_s) != NA_SCOMPLEX)
    rb_s = na_change_type(rb_s, NA_SCOMPLEX);
  s = NA_PTR_TYPE(rb_s, complex*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rb_x_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, complex*);
  MEMCPY(x_out__, x, complex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rb_y_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, complex*);
  MEMCPY(y_out__, y, complex, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rb_z_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  clar2v_(&n, x, y, z, &incx, c, s, &incc);

  return rb_ary_new3(3, rb_x, rb_y, rb_z);
}

void
init_lapack_clar2v(VALUE mLapack){
  rb_define_module_function(mLapack, "clar2v", rb_clar2v, -1);
}
