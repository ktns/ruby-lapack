#include "rb_lapack.h"

extern logical disnan_(doublereal *din);

static VALUE sHelp, sUsage;

static VALUE
rblapack_disnan(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_din;
  doublereal din; 
  VALUE rblapack___out__;
  logical __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.disnan( din, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      LOGICAL FUNCTION DISNAN( DIN )\n\n*  Purpose\n*  =======\n*\n*  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.\n*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the\n*  future.\n*\n\n*  Arguments\n*  =========\n*\n*  DIN     (input) DOUBLE PRECISION\n*          Input to test for NaN.\n*\n\n*  =====================================================================\n*\n*  .. External Functions ..\n      LOGICAL DLAISNAN\n      EXTERNAL DLAISNAN\n*  ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.disnan( din, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rblapack_din = argv[0];
  if (rb_options != Qnil) {
  }

  din = NUM2DBL(rblapack_din);

  __out__ = disnan_(&din);

  rblapack___out__ = __out__ ? Qtrue : Qfalse;
  return rblapack___out__;
}

void
init_lapack_disnan(VALUE mLapack){
  rb_define_module_function(mLapack, "disnan", rblapack_disnan, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
