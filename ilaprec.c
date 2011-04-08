#include "rb_lapack.h"

extern integer ilaprec_(char *prec);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ilaprec(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_prec;
  char prec; 
  VALUE rblapack___out__;
  integer __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilaprec( prec, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION ILAPREC( PREC )\n\n*  Purpose\n*  =======\n*\n*  This subroutine translated from a character string specifying an\n*  intermediate precision to the relevant BLAST-specified integer\n*  constant.\n*\n*  ILAPREC returns an INTEGER.  If ILAPREC < 0, then the input is not a\n*  character indicating a supported intermediate precision.  Otherwise\n*  ILAPREC returns the constant value corresponding to PREC.\n*\n\n*  Arguments\n*  =========\n*  PREC    (input) CHARACTER\n*          Specifies the form of the system of equations:\n*          = 'S':  Single\n*          = 'D':  Double\n*          = 'I':  Indigenous\n*          = 'X', 'E':  Extra\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilaprec( prec, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rblapack_prec = argv[0];
  if (rb_options != Qnil) {
  }

  prec = StringValueCStr(rblapack_prec)[0];

  __out__ = ilaprec_(&prec);

  rblapack___out__ = INT2NUM(__out__);
  return rblapack___out__;
}

void
init_lapack_ilaprec(VALUE mLapack){
  rb_define_module_function(mLapack, "ilaprec", rblapack_ilaprec, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
