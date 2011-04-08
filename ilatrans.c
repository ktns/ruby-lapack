#include "rb_lapack.h"

extern integer ilatrans_(char *trans);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ilatrans(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack___out__;
  integer __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilatrans( trans, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION ILATRANS( TRANS )\n\n*  Purpose\n*  =======\n*\n*  This subroutine translates from a character string specifying a\n*  transposition operation to the relevant BLAST-specified integer\n*  constant.\n*\n*  ILATRANS returns an INTEGER.  If ILATRANS < 0, then the input is not\n*  a character indicating a transposition operator.  Otherwise ILATRANS\n*  returns the constant value corresponding to TRANS.\n*\n\n*  Arguments\n*  =========\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  No transpose\n*          = 'T':  Transpose\n*          = 'C':  Conjugate transpose\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilatrans( trans, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rblapack_trans = argv[0];
  if (rb_options != Qnil) {
  }

  trans = StringValueCStr(rblapack_trans)[0];

  __out__ = ilatrans_(&trans);

  rblapack___out__ = INT2NUM(__out__);
  return rblapack___out__;
}

void
init_lapack_ilatrans(VALUE mLapack){
  rb_define_module_function(mLapack, "ilatrans", rblapack_ilatrans, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
