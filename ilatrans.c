#include "rb_lapack.h"

static VALUE
rb_ilatrans(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilatrans( trans)\n    or\n  NumRu::Lapack.ilatrans  # print help\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION ILATRANS( TRANS )\n\n*  Purpose\n*  =======\n*\n*  This subroutine translates from a character string specifying a\n*  transposition operation to the relevant BLAST-specified integer\n*  constant.\n*\n*  ILATRANS returns an INTEGER.  If ILATRANS < 0, then the input is not\n*  a character indicating a transposition operator.  Otherwise ILATRANS\n*  returns the constant value corresponding to TRANS.\n*\n\n*  Arguments\n*  =========\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  No transpose\n*          = 'T':  Transpose\n*          = 'C':  Conjugate transpose\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_trans = argv[0];

  trans = StringValueCStr(rb_trans)[0];

  __out__ = ilatrans_(&trans);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_ilatrans(VALUE mLapack){
  rb_define_module_function(mLapack, "ilatrans", rb_ilatrans, -1);
}
