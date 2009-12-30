#include "rb_lapack.h"

static VALUE
rb_slaruv(int argc, VALUE *argv, VALUE self){
  VALUE rb_iseed;
  integer *iseed; 
  VALUE rb_n;
  integer n; 
  VALUE rb_x;
  real *x; 
  VALUE rb_iseed_out__;
  integer *iseed_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, iseed = NumRu::Lapack.slaruv( iseed, n)\n    or\n  NumRu::Lapack.slaruv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARUV( ISEED, N, X )\n\n*  Purpose\n*  =======\n*\n*  SLARUV returns a vector of n random real numbers from a uniform (0,1)\n*  distribution (n <= 128).\n*\n*  This is an auxiliary routine called by SLARNV and CLARNV.\n*\n\n*  Arguments\n*  =========\n*\n*  ISEED   (input/output) INTEGER array, dimension (4)\n*          On entry, the seed of the random number generator; the array\n*          elements must be between 0 and 4095, and ISEED(4) must be\n*          odd.\n*          On exit, the seed is updated.\n*\n*  N       (input) INTEGER\n*          The number of random numbers to be generated. N <= 128.\n*\n*  X       (output) REAL array, dimension (N)\n*          The generated random numbers.\n*\n\n*  Further Details\n*  ===============\n*\n*  This routine uses a multiplicative congruential method with modulus\n*  2**48 and multiplier 33952834046453 (see G.S.Fishman,\n*  'Multiplicative congruential random number generators with modulus\n*  2**b: an exhaustive analysis for b = 32 and a partial analysis for\n*  b = 48', Math. Comp. 189, pp 331-344, 1990).\n*\n*  48-bit integers are stored in 4 integer array elements with 12 bits\n*  per element. Hence the routine is portable across machines with\n*  integers of 32 bits or more.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_iseed = argv[0];
  rb_n = argv[1];

  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_iseed))
    rb_raise(rb_eArgError, "iseed (1th argument) must be NArray");
  if (NA_RANK(rb_iseed) != 1)
    rb_raise(rb_eArgError, "rank of iseed (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iseed) != (4))
    rb_raise(rb_eRuntimeError, "shape 0 of iseed must be %d", 4);
  if (NA_TYPE(rb_iseed) != NA_LINT)
    rb_iseed = na_change_type(rb_iseed, NA_LINT);
  iseed = NA_PTR_TYPE(rb_iseed, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rb_x = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = DIM_LEN(4);
    rb_iseed_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iseed_out__ = NA_PTR_TYPE(rb_iseed_out__, integer*);
  MEMCPY(iseed_out__, iseed, integer, NA_TOTAL(rb_iseed));
  rb_iseed = rb_iseed_out__;
  iseed = iseed_out__;

  slaruv_(iseed, &n, x);

  return rb_ary_new3(2, rb_x, rb_iseed);
}

void
init_lapack_slaruv(VALUE mLapack){
  rb_define_module_function(mLapack, "slaruv", rb_slaruv, -1);
}
