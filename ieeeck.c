#include "rb_lapack.h"

extern integer ieeeck_(integer *ispec, real *zero, real *one);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ieeeck(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_ispec;
  integer ispec; 
  VALUE rblapack_zero;
  real zero; 
  VALUE rblapack_one;
  real one; 
  VALUE rblapack___out__;
  integer __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ieeeck( ispec, zero, one, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )\n\n*  Purpose\n*  =======\n*\n*  IEEECK is called from the ILAENV to verify that Infinity and\n*  possibly NaN arithmetic is safe (i.e. will not trap).\n*\n\n*  Arguments\n*  =========\n*\n*  ISPEC   (input) INTEGER\n*          Specifies whether to test just for inifinity arithmetic\n*          or whether to test for infinity and NaN arithmetic.\n*          = 0: Verify infinity arithmetic only.\n*          = 1: Verify infinity and NaN arithmetic.\n*\n*  ZERO    (input) REAL\n*          Must contain the value 0.0\n*          This is passed to prevent the compiler from optimizing\n*          away this code.\n*\n*  ONE     (input) REAL\n*          Must contain the value 1.0\n*          This is passed to prevent the compiler from optimizing\n*          away this code.\n*\n*  RETURN VALUE:  INTEGER\n*          = 0:  Arithmetic failed to produce the correct answers\n*          = 1:  Arithmetic produced the correct answers\n*\n*     .. Local Scalars ..\n      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,\n     $                   NEGZRO, NEWZRO, POSINF\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ieeeck( ispec, zero, one, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_ispec = argv[0];
  rblapack_zero = argv[1];
  rblapack_one = argv[2];
  if (rb_options != Qnil) {
  }

  one = (real)NUM2DBL(rblapack_one);
  ispec = NUM2INT(rblapack_ispec);
  zero = (real)NUM2DBL(rblapack_zero);

  __out__ = ieeeck_(&ispec, &zero, &one);

  rblapack___out__ = INT2NUM(__out__);
  return rblapack___out__;
}

void
init_lapack_ieeeck(VALUE mLapack){
  rb_define_module_function(mLapack, "ieeeck", rblapack_ieeeck, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
