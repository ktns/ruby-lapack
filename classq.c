#include "rb_lapack.h"

extern VOID classq_(integer *n, complex *x, integer *incx, real *scale, real *sumsq);

static VALUE sHelp, sUsage;

static VALUE
rblapack_classq(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  complex *x; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_scale;
  real scale; 
  VALUE rblapack_sumsq;
  real sumsq; 

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, sumsq = NumRu::Lapack.classq( x, incx, scale, sumsq, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLASSQ( N, X, INCX, SCALE, SUMSQ )\n\n*  Purpose\n*  =======\n*\n*  CLASSQ returns the values scl and ssq such that\n*\n*     ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,\n*\n*  where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is\n*  assumed to be at least unity and the value of ssq will then satisfy\n*\n*     1.0 .le. ssq .le. ( sumsq + 2*n ).\n*\n*  scale is assumed to be non-negative and scl returns the value\n*\n*     scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),\n*            i\n*\n*  scale and sumsq must be supplied in SCALE and SUMSQ respectively.\n*  SCALE and SUMSQ are overwritten by scl and ssq respectively.\n*\n*  The routine makes only one pass through the vector X.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of elements to be used from the vector X.\n*\n*  X       (input) COMPLEX array, dimension (N)\n*          The vector x as described above.\n*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of the vector X.\n*          INCX > 0.\n*\n*  SCALE   (input/output) REAL\n*          On entry, the value  scale  in the equation above.\n*          On exit, SCALE is overwritten with the value  scl .\n*\n*  SUMSQ   (input/output) REAL\n*          On entry, the value  sumsq  in the equation above.\n*          On exit, SUMSQ is overwritten with the value  ssq .\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, sumsq = NumRu::Lapack.classq( x, incx, scale, sumsq, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_x = argv[0];
  rblapack_incx = argv[1];
  rblapack_scale = argv[2];
  rblapack_sumsq = argv[3];
  if (rb_options != Qnil) {
  }

  scale = (real)NUM2DBL(rblapack_scale);
  sumsq = (real)NUM2DBL(rblapack_sumsq);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_SCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, complex*);
  incx = NUM2INT(rblapack_incx);

  classq_(&n, x, &incx, &scale, &sumsq);

  rblapack_scale = rb_float_new((double)scale);
  rblapack_sumsq = rb_float_new((double)sumsq);
  return rb_ary_new3(2, rblapack_scale, rblapack_sumsq);
}

void
init_lapack_classq(VALUE mLapack){
  rb_define_module_function(mLapack, "classq", rblapack_classq, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
