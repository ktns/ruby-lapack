#include "rb_lapack.h"

extern integer ilauplo_(char *uplo);

static VALUE
rb_ilauplo(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilauplo( uplo)\n    or\n  NumRu::Lapack.ilauplo  # print help\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION ILAUPLO( UPLO )\n\n*  Purpose\n*  =======\n*\n*  This subroutine translated from a character string specifying a\n*  upper- or lower-triangular matrix to the relevant BLAST-specified\n*  integer constant.\n*\n*  ILAUPLO returns an INTEGER.  If ILAUPLO < 0, then the input is not\n*  a character indicating an upper- or lower-triangular matrix.\n*  Otherwise ILAUPLO returns the constant value corresponding to UPLO.\n*\n\n*  Arguments\n*  =========\n*  UPLO    (input) CHARACTER\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_uplo = argv[0];

  uplo = StringValueCStr(rb_uplo)[0];

  __out__ = ilauplo_(&uplo);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_ilauplo(VALUE mLapack){
  rb_define_module_function(mLapack, "ilauplo", rb_ilauplo, -1);
}
