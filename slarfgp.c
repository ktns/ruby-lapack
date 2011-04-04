#include "rb_lapack.h"

extern VOID slarfgp_(integer *n, real *alpha, real *x, integer *incx, real *tau);

static VALUE
rb_slarfgp(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_alpha;
  real alpha; 
  VALUE rb_x;
  real *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_tau;
  real tau; 
  VALUE rb_x_out__;
  real *x_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, alpha, x = NumRu::Lapack.slarfgp( n, alpha, x, incx)\n    or\n  NumRu::Lapack.slarfgp  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_n = argv[0];
  rb_alpha = argv[1];
  rb_x = argv[2];
  rb_incx = argv[3];

  alpha = (real)NUM2DBL(rb_alpha);
  n = NUM2INT(rb_n);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (3th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-2)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-2)*abs(incx));
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = 1+(n-2)*abs(incx);
    rb_x_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;

  slarfgp_(&n, &alpha, x, &incx, &tau);

  rb_tau = rb_float_new((double)tau);
  rb_alpha = rb_float_new((double)alpha);
  return rb_ary_new3(3, rb_tau, rb_alpha, rb_x);
}

void
init_lapack_slarfgp(VALUE mLapack){
  rb_define_module_function(mLapack, "slarfgp", rb_slarfgp, -1);
}
