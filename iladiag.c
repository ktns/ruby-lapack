#include "rb_lapack.h"

static VALUE
rb_iladiag(int argc, VALUE *argv, VALUE self){
  VALUE rb_diag;
  char diag; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.iladiag( diag)\n    or\n  NumRu::Lapack.iladiag  # print help\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION ILADIAG( DIAG )\n\n*  Purpose\n*  =======\n*\n*  This subroutine translated from a character string specifying if a\n*  matrix has unit diagonal or not to the relevant BLAST-specified\n*  integer constant.\n*\n*  ILADIAG returns an INTEGER.  If ILADIAG < 0, then the input is not a\n*  character indicating a unit or non-unit diagonal.  Otherwise ILADIAG\n*  returns the constant value corresponding to DIAG.\n*\n\n*  Arguments\n*  =========\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_diag = argv[0];

  diag = StringValueCStr(rb_diag)[0];

  __out__ = iladiag_(&diag);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_iladiag(VALUE mLapack){
  rb_define_module_function(mLapack, "iladiag", rb_iladiag, -1);
}
