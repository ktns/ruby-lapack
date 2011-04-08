#include "rb_lapack.h"

extern VOID slaruv_(integer *iseed, integer *n, real *x);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaruv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_iseed;
  integer *iseed; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_x;
  real *x; 
  VALUE rblapack_iseed_out__;
  integer *iseed_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, iseed = NumRu::Lapack.slaruv( iseed, n, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARUV( ISEED, N, X )\n\n*  Purpose\n*  =======\n*\n*  SLARUV returns a vector of n random real numbers from a uniform (0,1)\n*  distribution (n <= 128).\n*\n*  This is an auxiliary routine called by SLARNV and CLARNV.\n*\n\n*  Arguments\n*  =========\n*\n*  ISEED   (input/output) INTEGER array, dimension (4)\n*          On entry, the seed of the random number generator; the array\n*          elements must be between 0 and 4095, and ISEED(4) must be\n*          odd.\n*          On exit, the seed is updated.\n*\n*  N       (input) INTEGER\n*          The number of random numbers to be generated. N <= 128.\n*\n*  X       (output) REAL array, dimension (N)\n*          The generated random numbers.\n*\n\n*  Further Details\n*  ===============\n*\n*  This routine uses a multiplicative congruential method with modulus\n*  2**48 and multiplier 33952834046453 (see G.S.Fishman,\n*  'Multiplicative congruential random number generators with modulus\n*  2**b: an exhaustive analysis for b = 32 and a partial analysis for\n*  b = 48', Math. Comp. 189, pp 331-344, 1990).\n*\n*  48-bit integers are stored in 4 integer array elements with 12 bits\n*  per element. Hence the routine is portable across machines with\n*  integers of 32 bits or more.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, iseed = NumRu::Lapack.slaruv( iseed, n, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_iseed = argv[0];
  rblapack_n = argv[1];
  if (rb_options != Qnil) {
  }

  n = NUM2INT(rblapack_n);
  if (!NA_IsNArray(rblapack_iseed))
    rb_raise(rb_eArgError, "iseed (1th argument) must be NArray");
  if (NA_RANK(rblapack_iseed) != 1)
    rb_raise(rb_eArgError, "rank of iseed (1th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_iseed) != (4))
    rb_raise(rb_eRuntimeError, "shape 0 of iseed must be %d", 4);
  if (NA_TYPE(rblapack_iseed) != NA_LINT)
    rblapack_iseed = na_change_type(rblapack_iseed, NA_LINT);
  iseed = NA_PTR_TYPE(rblapack_iseed, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rblapack_x = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, real*);
  {
    int shape[1];
    shape[0] = 4;
    rblapack_iseed_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iseed_out__ = NA_PTR_TYPE(rblapack_iseed_out__, integer*);
  MEMCPY(iseed_out__, iseed, integer, NA_TOTAL(rblapack_iseed));
  rblapack_iseed = rblapack_iseed_out__;
  iseed = iseed_out__;

  slaruv_(iseed, &n, x);

  return rb_ary_new3(2, rblapack_x, rblapack_iseed);
}

void
init_lapack_slaruv(VALUE mLapack){
  rb_define_module_function(mLapack, "slaruv", rblapack_slaruv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
