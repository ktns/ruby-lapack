#include "rb_lapack.h"

extern VOID dlarfgp_(integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *tau);

static VALUE
rb_dlarfgp(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_alpha;
  doublereal alpha; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_tau;
  doublereal tau; 
  VALUE rb_x_out__;
  doublereal *x_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, alpha, x = NumRu::Lapack.dlarfgp( n, alpha, x, incx)\n    or\n  NumRu::Lapack.dlarfgp  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARFGP( N, ALPHA, X, INCX, TAU )\n\n*  Purpose\n*  =======\n*\n*  DLARFGP generates a real elementary reflector H of order n, such\n*  that\n*\n*        H * ( alpha ) = ( beta ),   H' * H = I.\n*            (   x   )   (   0  )\n*\n*  where alpha and beta are scalars, beta is non-negative, and x is\n*  an (n-1)-element real vector.  H is represented in the form\n*\n*        H = I - tau * ( 1 ) * ( 1 v' ) ,\n*                      ( v )\n*\n*  where tau is a real scalar and v is a real (n-1)-element\n*  vector.\n*\n*  If the elements of x are all zero, then tau = 0 and H is taken to be\n*  the unit matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the elementary reflector.\n*\n*  ALPHA   (input/output) DOUBLE PRECISION\n*          On entry, the value alpha.\n*          On exit, it is overwritten with the value beta.\n*\n*  X       (input/output) DOUBLE PRECISION array, dimension\n*                         (1+(N-2)*abs(INCX))\n*          On entry, the vector x.\n*          On exit, it is overwritten with the vector v.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X. INCX > 0.\n*\n*  TAU     (output) DOUBLE PRECISION\n*          The value tau.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_n = argv[0];
  rb_alpha = argv[1];
  rb_x = argv[2];
  rb_incx = argv[3];

  alpha = NUM2DBL(rb_alpha);
  n = NUM2INT(rb_n);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (3th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-2)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-2)*abs(incx));
  if (NA_TYPE(rb_x) != NA_DFLOAT)
    rb_x = na_change_type(rb_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = 1+(n-2)*abs(incx);
    rb_x_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublereal*);
  MEMCPY(x_out__, x, doublereal, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;

  dlarfgp_(&n, &alpha, x, &incx, &tau);

  rb_tau = rb_float_new((double)tau);
  rb_alpha = rb_float_new((double)alpha);
  return rb_ary_new3(3, rb_tau, rb_alpha, rb_x);
}

void
init_lapack_dlarfgp(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarfgp", rb_dlarfgp, -1);
}
