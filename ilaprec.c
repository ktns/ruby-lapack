#include "rb_lapack.h"

static VALUE
rb_ilaprec(int argc, VALUE *argv, VALUE self){
  VALUE rb_prec;
  char prec; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilaprec( prec)\n    or\n  NumRu::Lapack.ilaprec  # print help\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION ILAPREC( PREC )\n\n*  Purpose\n*  =======\n*\n*  This subroutine translated from a character string specifying an\n*  intermediate precision to the relevant BLAST-specified integer\n*  constant.\n*\n*  ILAPREC returns an INTEGER.  If ILAPREC < 0, then the input is not a\n*  character indicating a supported intermediate precision.  Otherwise\n*  ILAPREC returns the constant value corresponding to PREC.\n*\n\n*  Arguments\n*  =========\n*  PREC   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'S':  Single\n*          = 'D':  Double\n*          = 'I':  Indigenous\n*          = 'X', 'E':  Extra\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_prec = argv[0];

  prec = StringValueCStr(rb_prec)[0];

  __out__ = ilaprec_(&prec);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_ilaprec(VALUE mLapack){
  rb_define_module_function(mLapack, "ilaprec", rb_ilaprec, -1);
}
