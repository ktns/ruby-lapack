#include "rb_lapack.h"

extern VOID slabad_(real *small, real *large);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slabad(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_small;
  real small; 
  VALUE rblapack_large;
  real large; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  small, large = NumRu::Lapack.slabad( small, large, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLABAD( SMALL, LARGE )\n\n*  Purpose\n*  =======\n*\n*  SLABAD takes as input the values computed by SLAMCH for underflow and\n*  overflow, and returns the square root of each of these values if the\n*  log of LARGE is sufficiently large.  This subroutine is intended to\n*  identify machines with a large exponent range, such as the Crays, and\n*  redefine the underflow and overflow limits to be the square roots of\n*  the values computed by SLAMCH.  This subroutine is needed because\n*  SLAMCH does not compensate for poor arithmetic in the upper half of\n*  the exponent range, as is found on a Cray.\n*\n\n*  Arguments\n*  =========\n*\n*  SMALL   (input/output) REAL\n*          On entry, the underflow threshold as computed by SLAMCH.\n*          On exit, if LOG10(LARGE) is sufficiently large, the square\n*          root of SMALL, otherwise unchanged.\n*\n*  LARGE   (input/output) REAL\n*          On entry, the overflow threshold as computed by SLAMCH.\n*          On exit, if LOG10(LARGE) is sufficiently large, the square\n*          root of LARGE, otherwise unchanged.\n*\n\n*  =====================================================================\n*\n*     .. Intrinsic Functions ..\n      INTRINSIC          LOG10, SQRT\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  small, large = NumRu::Lapack.slabad( small, large, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_small = argv[0];
  rblapack_large = argv[1];
  if (rb_options != Qnil) {
  }

  large = (real)NUM2DBL(rblapack_large);
  small = (real)NUM2DBL(rblapack_small);

  slabad_(&small, &large);

  rblapack_small = rb_float_new((double)small);
  rblapack_large = rb_float_new((double)large);
  return rb_ary_new3(2, rblapack_small, rblapack_large);
}

void
init_lapack_slabad(VALUE mLapack){
  rb_define_module_function(mLapack, "slabad", rblapack_slabad, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
