#include "rb_lapack.h"

static VALUE
rb_clarfg(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_alpha;
  complex alpha; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_tau;
  complex tau; 
  VALUE rb_x_out__;
  complex *x_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, alpha, x = NumRu::Lapack.clarfg( n, alpha, x, incx)\n    or\n  NumRu::Lapack.clarfg  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU )\n\n*  Purpose\n*  =======\n*\n*  CLARFG generates a complex elementary reflector H of order n, such\n*  that\n*\n*        H' * ( alpha ) = ( beta ),   H' * H = I.\n*             (   x   )   (   0  )\n*\n*  where alpha and beta are scalars, with beta real, and x is an\n*  (n-1)-element complex vector. H is represented in the form\n*\n*        H = I - tau * ( 1 ) * ( 1 v' ) ,\n*                      ( v )\n*\n*  where tau is a complex scalar and v is a complex (n-1)-element\n*  vector. Note that H is not hermitian.\n*\n*  If the elements of x are all zero and alpha is real, then tau = 0\n*  and H is taken to be the unit matrix.\n*\n*  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the elementary reflector.\n*\n*  ALPHA   (input/output) COMPLEX\n*          On entry, the value alpha.\n*          On exit, it is overwritten with the value beta.\n*\n*  X       (input/output) COMPLEX array, dimension\n*                         (1+(N-2)*abs(INCX))\n*          On entry, the vector x.\n*          On exit, it is overwritten with the vector v.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X. INCX > 0.\n*\n*  TAU     (output) COMPLEX\n*          The value tau.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_n = argv[0];
  rb_alpha = argv[1];
  rb_x = argv[2];
  rb_incx = argv[3];

  n = NUM2INT(rb_n);
  alpha.r = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (3th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-2)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-2)*abs(incx));
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  {
    int shape[1];
    shape[0] = 1+(n-2)*abs(incx);
    rb_x_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, complex*);
  MEMCPY(x_out__, x, complex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;

  clarfg_(&n, &alpha, x, &incx, &tau);

  rb_tau = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(tau.r)), rb_float_new((double)(tau.i)));
  rb_alpha = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(alpha.r)), rb_float_new((double)(alpha.i)));
  return rb_ary_new3(3, rb_tau, rb_alpha, rb_x);
}

void
init_lapack_clarfg(VALUE mLapack){
  rb_define_module_function(mLapack, "clarfg", rb_clarfg, -1);
}
