#include "rb_lapack.h"

extern VOID slassq_(integer *n, real *x, integer *incx, real *scale, real *sumsq);

static VALUE
rb_slassq(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  real *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_scale;
  real scale; 
  VALUE rb_sumsq;
  real sumsq; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, sumsq = NumRu::Lapack.slassq( x, incx, scale, sumsq)\n    or\n  NumRu::Lapack.slassq  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ )\n\n*  Purpose\n*  =======\n*\n*  SLASSQ  returns the values  scl  and  smsq  such that\n*\n*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,\n*\n*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is\n*  assumed to be non-negative and  scl  returns the value\n*\n*     scl = max( scale, abs( x( i ) ) ).\n*\n*  scale and sumsq must be supplied in SCALE and SUMSQ and\n*  scl and smsq are overwritten on SCALE and SUMSQ respectively.\n*\n*  The routine makes only one pass through the vector x.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of elements to be used from the vector X.\n*\n*  X       (input) REAL array, dimension (N)\n*          The vector for which a scaled sum of squares is computed.\n*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of the vector X.\n*          INCX > 0.\n*\n*  SCALE   (input/output) REAL\n*          On entry, the value  scale  in the equation above.\n*          On exit, SCALE is overwritten with  scl , the scaling factor\n*          for the sum of squares.\n*\n*  SUMSQ   (input/output) REAL\n*          On entry, the value  sumsq  in the equation above.\n*          On exit, SUMSQ is overwritten with  smsq , the basic sum of\n*          squares from which  scl  has been factored out.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_x = argv[0];
  rb_incx = argv[1];
  rb_scale = argv[2];
  rb_sumsq = argv[3];

  scale = (real)NUM2DBL(rb_scale);
  sumsq = (real)NUM2DBL(rb_sumsq);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  incx = NUM2INT(rb_incx);

  slassq_(&n, x, &incx, &scale, &sumsq);

  rb_scale = rb_float_new((double)scale);
  rb_sumsq = rb_float_new((double)sumsq);
  return rb_ary_new3(2, rb_scale, rb_sumsq);
}

void
init_lapack_slassq(VALUE mLapack){
  rb_define_module_function(mLapack, "slassq", rb_slassq, -1);
}
