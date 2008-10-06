#include "rb_lapack.h"

static VALUE
rb_slabad(int argc, VALUE *argv, VALUE self){
  VALUE rb_small;
  real small; 
  VALUE rb_large;
  real large; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  small, large = NumRu::Lapack.slabad( small, large)\n    or\n  NumRu::Lapack.slabad  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLABAD( SMALL, LARGE )\n\n*  Purpose\n*  =======\n*\n*  SLABAD takes as input the values computed by SLAMCH for underflow and\n*  overflow, and returns the square root of each of these values if the\n*  log of LARGE is sufficiently large.  This subroutine is intended to\n*  identify machines with a large exponent range, such as the Crays, and\n*  redefine the underflow and overflow limits to be the square roots of\n*  the values computed by SLAMCH.  This subroutine is needed because\n*  SLAMCH does not compensate for poor arithmetic in the upper half of\n*  the exponent range, as is found on a Cray.\n*\n\n*  Arguments\n*  =========\n*\n*  SMALL   (input/output) REAL\n*          On entry, the underflow threshold as computed by SLAMCH.\n*          On exit, if LOG10(LARGE) is sufficiently large, the square\n*          root of SMALL, otherwise unchanged.\n*\n*  LARGE   (input/output) REAL\n*          On entry, the overflow threshold as computed by SLAMCH.\n*          On exit, if LOG10(LARGE) is sufficiently large, the square\n*          root of LARGE, otherwise unchanged.\n*\n\n*  =====================================================================\n*\n*     .. Intrinsic Functions ..\n      INTRINSIC          LOG10, SQRT\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_small = argv[0];
  rb_large = argv[1];

  small = (real)NUM2DBL(rb_small);
  large = (real)NUM2DBL(rb_large);

  slabad_(&small, &large);

  rb_small = rb_float_new((double)small);
  rb_large = rb_float_new((double)large);
  return rb_ary_new3(2, rb_small, rb_large);
}

void
init_lapack_slabad(VALUE mLapack){
  rb_define_module_function(mLapack, "slabad", rb_slabad, -1);
}
