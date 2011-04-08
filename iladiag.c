#include "rb_lapack.h"

extern integer iladiag_(char *diag);

static VALUE sHelp, sUsage;

static VALUE
rblapack_iladiag(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_diag;
  char diag; 
  VALUE rblapack___out__;
  integer __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.iladiag( diag, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION ILADIAG( DIAG )\n\n*  Purpose\n*  =======\n*\n*  This subroutine translated from a character string specifying if a\n*  matrix has unit diagonal or not to the relevant BLAST-specified\n*  integer constant.\n*\n*  ILADIAG returns an INTEGER.  If ILADIAG < 0, then the input is not a\n*  character indicating a unit or non-unit diagonal.  Otherwise ILADIAG\n*  returns the constant value corresponding to DIAG.\n*\n\n*  Arguments\n*  =========\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.iladiag( diag, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rblapack_diag = argv[0];
  if (rb_options != Qnil) {
  }

  diag = StringValueCStr(rblapack_diag)[0];

  __out__ = iladiag_(&diag);

  rblapack___out__ = INT2NUM(__out__);
  return rblapack___out__;
}

void
init_lapack_iladiag(VALUE mLapack){
  rb_define_module_function(mLapack, "iladiag", rblapack_iladiag, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
