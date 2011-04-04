#include "rb_lapack.h"

extern integer ieeeck_(integer *ispec, real *zero, real *one);

static VALUE
rb_ieeeck(int argc, VALUE *argv, VALUE self){
  VALUE rb_ispec;
  integer ispec; 
  VALUE rb_zero;
  real zero; 
  VALUE rb_one;
  real one; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ieeeck( ispec, zero, one)\n    or\n  NumRu::Lapack.ieeeck  # print help\n\n\nFORTRAN MANUAL\n      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )\n\n*  Purpose\n*  =======\n*\n*  IEEECK is called from the ILAENV to verify that Infinity and\n*  possibly NaN arithmetic is safe (i.e. will not trap).\n*\n\n*  Arguments\n*  =========\n*\n*  ISPEC   (input) INTEGER\n*          Specifies whether to test just for inifinity arithmetic\n*          or whether to test for infinity and NaN arithmetic.\n*          = 0: Verify infinity arithmetic only.\n*          = 1: Verify infinity and NaN arithmetic.\n*\n*  ZERO    (input) REAL\n*          Must contain the value 0.0\n*          This is passed to prevent the compiler from optimizing\n*          away this code.\n*\n*  ONE     (input) REAL\n*          Must contain the value 1.0\n*          This is passed to prevent the compiler from optimizing\n*          away this code.\n*\n*  RETURN VALUE:  INTEGER\n*          = 0:  Arithmetic failed to produce the correct answers\n*          = 1:  Arithmetic produced the correct answers\n*\n*     .. Local Scalars ..\n      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,\n     $                   NEGZRO, NEWZRO, POSINF\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_ispec = argv[0];
  rb_zero = argv[1];
  rb_one = argv[2];

  one = (real)NUM2DBL(rb_one);
  ispec = NUM2INT(rb_ispec);
  zero = (real)NUM2DBL(rb_zero);

  __out__ = ieeeck_(&ispec, &zero, &one);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_ieeeck(VALUE mLapack){
  rb_define_module_function(mLapack, "ieeeck", rb_ieeeck, -1);
}
