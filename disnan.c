#include "rb_lapack.h"

static VALUE
rb_disnan(int argc, VALUE *argv, VALUE self){
  VALUE rb_din;
  doublereal din; 
  VALUE rb___out__;
  logical __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.disnan( din)\n    or\n  NumRu::Lapack.disnan  # print help\n\n\nFORTRAN MANUAL\n      LOGICAL FUNCTION DISNAN(DIN)\n\n*  Purpose\n*  =======\n*\n*  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.\n*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the\n*  future.\n*\n\n*  Arguments\n*  =========\n*\n*  DIN      (input) DOUBLE PRECISION\n*          Input to test for NaN.\n*\n\n*  =====================================================================\n*\n*  .. External Functions ..\n      LOGICAL DLAISNAN\n      EXTERNAL DLAISNAN\n*  ..\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_din = argv[0];

  din = NUM2DBL(rb_din);

  __out__ = disnan_(&din);

  rb___out__ = __out__ ? Qtrue : Qfalse;
  return rb___out__;
}

void
init_lapack_disnan(VALUE mLapack){
  rb_define_module_function(mLapack, "disnan", rb_disnan, -1);
}
