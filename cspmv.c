#include "rb_lapack.h"

extern VOID cspmv_(char *uplo, integer *n, complex *alpha, complex *ap, complex *x, integer *incx, complex *beta, complex *y, integer *incy);

static VALUE
rb_cspmv(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_alpha;
  complex alpha; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_beta;
  complex beta; 
  VALUE rb_y;
  complex *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_y_out__;
  complex *y_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  y = NumRu::Lapack.cspmv( uplo, n, alpha, ap, x, incx, beta, y, incy)\n    or\n  NumRu::Lapack.cspmv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_uplo = argv[0];
  rb_n = argv[1];
  rb_alpha = argv[2];
  rb_ap = argv[3];
  rb_x = argv[4];
  rb_incx = argv[5];
  rb_beta = argv[6];
  rb_y = argv[7];
  rb_incy = argv[8];

  uplo = StringValueCStr(rb_uplo)[0];
  alpha.r = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  n = NUM2INT(rb_n);
  incy = NUM2INT(rb_incy);
  incx = NUM2INT(rb_incx);
  beta.r = (real)NUM2DBL(rb_funcall(rb_beta, rb_intern("real"), 0));
  beta.i = (real)NUM2DBL(rb_funcall(rb_beta, rb_intern("imag"), 0));
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (( n*( n + 1 ) )/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", ( n*( n + 1 ) )/2);
  if (NA_TYPE(rb_ap) != NA_SCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, complex*);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (8th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rb_y) != NA_SCOMPLEX)
    rb_y = na_change_type(rb_y, NA_SCOMPLEX);
  y = NA_PTR_TYPE(rb_y, complex*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (5th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1 + ( n - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1 + ( n - 1 )*abs( incx ));
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  {
    int shape[1];
    shape[0] = 1 + ( n - 1 )*abs( incy );
    rb_y_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, complex*);
  MEMCPY(y_out__, y, complex, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  cspmv_(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);

  return rb_y;
}

void
init_lapack_cspmv(VALUE mLapack){
  rb_define_module_function(mLapack, "cspmv", rb_cspmv, -1);
}
