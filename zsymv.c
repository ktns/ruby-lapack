#include "rb_lapack.h"

extern VOID zsymv_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);

static VALUE
rb_zsymv(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_alpha;
  doublecomplex alpha; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_beta;
  doublecomplex beta; 
  VALUE rb_y;
  doublecomplex *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_y_out__;
  doublecomplex *y_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  y = NumRu::Lapack.zsymv( uplo, alpha, a, x, incx, beta, y, incy)\n    or\n  NumRu::Lapack.zsymv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_uplo = argv[0];
  rb_alpha = argv[1];
  rb_a = argv[2];
  rb_x = argv[3];
  rb_incx = argv[4];
  rb_beta = argv[5];
  rb_y = argv[6];
  rb_incy = argv[7];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  uplo = StringValueCStr(rb_uplo)[0];
  alpha.r = NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  beta.r = NUM2DBL(rb_funcall(rb_beta, rb_intern("real"), 0));
  beta.i = NUM2DBL(rb_funcall(rb_beta, rb_intern("imag"), 0));
  incy = NUM2INT(rb_incy);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (7th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rb_y) != NA_DCOMPLEX)
    rb_y = na_change_type(rb_y, NA_DCOMPLEX);
  y = NA_PTR_TYPE(rb_y, doublecomplex*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (4th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1 + ( n - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1 + ( n - 1 )*abs( incx ));
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  {
    int shape[1];
    shape[0] = 1 + ( n - 1 )*abs( incy );
    rb_y_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, doublecomplex*);
  MEMCPY(y_out__, y, doublecomplex, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  zsymv_(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

  return rb_y;
}

void
init_lapack_zsymv(VALUE mLapack){
  rb_define_module_function(mLapack, "zsymv", rb_zsymv, -1);
}
