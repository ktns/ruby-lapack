#include "rb_lapack.h"

extern VOID chla_transtype_(char *__out__, integer *trans);
static VALUE
rb_chla_transtype(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  integer trans; 
  VALUE rb___out__;
  char __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.chla_transtype( trans)\n    or\n  NumRu::Lapack.chla_transtype  # print help\n\n\nFORTRAN MANUAL\n      CHARACTER*1 FUNCTION CHLA_TRANSTYPE( TRANS )\n\n*  Purpose\n*  =======\n*\n*  This subroutine translates from a BLAST-specified integer constant to\n*  the character string specifying a transposition operation.\n*\n*  CHLA_TRANSTYPE returns an CHARACTER*1.  If CHLA_TRANSTYPE is 'X',\n*  then input is not an integer indicating a transposition operator.\n*  Otherwise CHLA_TRANSTYPE returns the constant value corresponding to\n*  TRANS.\n*\n\n*  Arguments\n*  =========\n*  TRANS   (input) INTEGER\n*          Specifies the form of the system of equations:\n*          = BLAS_NO_TRANS   = 111 :  No Transpose\n*          = BLAS_TRANS      = 112 :  Transpose\n*          = BLAS_CONJ_TRANS = 113 :  Conjugate Transpose\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_trans = argv[0];

  trans = NUM2INT(rb_trans);

  chla_transtype_(&__out__, &trans);

  rb___out__ = rb_str_new(&__out__,1);
  return rb___out__;
}

void
init_lapack_chla_transtype(VALUE mLapack){
  rb_define_module_function(mLapack, "chla_transtype", rb_chla_transtype, -1);
}
