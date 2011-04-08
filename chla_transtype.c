#include "rb_lapack.h"

extern VOID chla_transtype_(char *__out__, integer *trans);

static VALUE sHelp, sUsage;

static VALUE
rblapack_chla_transtype(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  integer trans; 
  VALUE rblapack___out__;
  char __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.chla_transtype( trans, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      CHARACTER*1 FUNCTION CHLA_TRANSTYPE( TRANS )\n\n*  Purpose\n*  =======\n*\n*  This subroutine translates from a BLAST-specified integer constant to\n*  the character string specifying a transposition operation.\n*\n*  CHLA_TRANSTYPE returns an CHARACTER*1.  If CHLA_TRANSTYPE is 'X',\n*  then input is not an integer indicating a transposition operator.\n*  Otherwise CHLA_TRANSTYPE returns the constant value corresponding to\n*  TRANS.\n*\n\n*  Arguments\n*  =========\n*  TRANS   (input) INTEGER\n*          Specifies the form of the system of equations:\n*          = BLAS_NO_TRANS   = 111 :  No Transpose\n*          = BLAS_TRANS      = 112 :  Transpose\n*          = BLAS_CONJ_TRANS = 113 :  Conjugate Transpose\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.chla_transtype( trans, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rblapack_trans = argv[0];
  if (rb_options != Qnil) {
  }

  trans = NUM2INT(rblapack_trans);

  chla_transtype_(&__out__, &trans);

  rblapack___out__ = rb_str_new(&__out__,1);
  return rblapack___out__;
}

void
init_lapack_chla_transtype(VALUE mLapack){
  rb_define_module_function(mLapack, "chla_transtype", rblapack_chla_transtype, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
