#include "rb_lapack.h"

extern VOID dlarfg_(integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *tau);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlarfg(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_alpha;
  doublereal alpha; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_tau;
  doublereal tau; 
  VALUE rblapack_x_out__;
  doublereal *x_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, alpha, x = NumRu::Lapack.dlarfg( n, alpha, x, incx, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )\n\n*  Purpose\n*  =======\n*\n*  DLARFG generates a real elementary reflector H of order n, such\n*  that\n*\n*        H * ( alpha ) = ( beta ),   H' * H = I.\n*            (   x   )   (   0  )\n*\n*  where alpha and beta are scalars, and x is an (n-1)-element real\n*  vector. H is represented in the form\n*\n*        H = I - tau * ( 1 ) * ( 1 v' ) ,\n*                      ( v )\n*\n*  where tau is a real scalar and v is a real (n-1)-element\n*  vector.\n*\n*  If the elements of x are all zero, then tau = 0 and H is taken to be\n*  the unit matrix.\n*\n*  Otherwise  1 <= tau <= 2.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the elementary reflector.\n*\n*  ALPHA   (input/output) DOUBLE PRECISION\n*          On entry, the value alpha.\n*          On exit, it is overwritten with the value beta.\n*\n*  X       (input/output) DOUBLE PRECISION array, dimension\n*                         (1+(N-2)*abs(INCX))\n*          On entry, the vector x.\n*          On exit, it is overwritten with the vector v.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X. INCX > 0.\n*\n*  TAU     (output) DOUBLE PRECISION\n*          The value tau.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, alpha, x = NumRu::Lapack.dlarfg( n, alpha, x, incx, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_n = argv[0];
  rblapack_alpha = argv[1];
  rblapack_x = argv[2];
  rblapack_incx = argv[3];
  if (rb_options != Qnil) {
  }

  alpha = NUM2DBL(rblapack_alpha);
  n = NUM2INT(rblapack_n);
  incx = NUM2INT(rblapack_incx);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (3th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (1+(n-2)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-2)*abs(incx));
  if (NA_TYPE(rblapack_x) != NA_DFLOAT)
    rblapack_x = na_change_type(rblapack_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rblapack_x, doublereal*);
  {
    int shape[1];
    shape[0] = 1+(n-2)*abs(incx);
    rblapack_x_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, doublereal*);
  MEMCPY(x_out__, x, doublereal, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;

  dlarfg_(&n, &alpha, x, &incx, &tau);

  rblapack_tau = rb_float_new((double)tau);
  rblapack_alpha = rb_float_new((double)alpha);
  return rb_ary_new3(3, rblapack_tau, rblapack_alpha, rblapack_x);
}

void
init_lapack_dlarfg(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarfg", rblapack_dlarfg, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
