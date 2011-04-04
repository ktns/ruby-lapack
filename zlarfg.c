#include "rb_lapack.h"

extern VOID zlarfg_(integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *tau);

static VALUE
rb_zlarfg(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_alpha;
  doublecomplex alpha; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_tau;
  doublecomplex tau; 
  VALUE rb_x_out__;
  doublecomplex *x_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, alpha, x = NumRu::Lapack.zlarfg( n, alpha, x, incx)\n    or\n  NumRu::Lapack.zlarfg  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_n = argv[0];
  rb_alpha = argv[1];
  rb_x = argv[2];
  rb_incx = argv[3];

  alpha.r = NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  n = NUM2INT(rb_n);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (3th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-2)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-2)*abs(incx));
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  {
    int shape[1];
    shape[0] = 1+(n-2)*abs(incx);
    rb_x_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;

  zlarfg_(&n, &alpha, x, &incx, &tau);

  rb_tau = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(tau.r)), rb_float_new((double)(tau.i)));
  rb_alpha = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(alpha.r)), rb_float_new((double)(alpha.i)));
  return rb_ary_new3(3, rb_tau, rb_alpha, rb_x);
}

void
init_lapack_zlarfg(VALUE mLapack){
  rb_define_module_function(mLapack, "zlarfg", rb_zlarfg, -1);
}
