#include "rb_lapack.h"

extern logical sisnan_(real *sin);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sisnan(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_sin;
  real sin; 
  VALUE rblapack___out__;
  logical __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.sisnan( sin, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      LOGICAL FUNCTION SISNAN( SIN )\n\n*  Purpose\n*  =======\n*\n*  SISNAN returns .TRUE. if its argument is NaN, and .FALSE.\n*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the\n*  future.\n*\n\n*  Arguments\n*  =========\n*\n*  SIN     (input) REAL\n*          Input to test for NaN.\n*\n\n*  =====================================================================\n*\n*  .. External Functions ..\n      LOGICAL SLAISNAN\n      EXTERNAL SLAISNAN\n*  ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.sisnan( sin, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rblapack_sin = argv[0];
  if (rb_options != Qnil) {
  }

  sin = (real)NUM2DBL(rblapack_sin);

  __out__ = sisnan_(&sin);

  rblapack___out__ = __out__ ? Qtrue : Qfalse;
  return rblapack___out__;
}

void
init_lapack_sisnan(VALUE mLapack){
  rb_define_module_function(mLapack, "sisnan", rblapack_sisnan, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
