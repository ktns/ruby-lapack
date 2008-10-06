#include "rb_lapack.h"

extern VOID sisnan_(logical *__out__, real *sin);
static VALUE
rb_sisnan(int argc, VALUE *argv, VALUE self){
  VALUE rb_sin;
  real sin; 
  VALUE rb___out__;
  logical __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.sisnan( sin)\n    or\n  NumRu::Lapack.sisnan  # print help\n\n\nFORTRAN MANUAL\n      LOGICAL FUNCTION SISNAN(SIN)\n\n*  Purpose\n*  =======\n*\n*  SISNAN returns .TRUE. if its argument is NaN, and .FALSE.\n*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the\n*  future.\n*\n\n*  Arguments\n*  =========\n*\n*  SIN      (input) REAL\n*          Input to test for NaN.\n*\n\n*  =====================================================================\n*\n*  .. External Functions ..\n      LOGICAL SLAISNAN\n      EXTERNAL SLAISNAN\n*  ..\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_sin = argv[0];

  sin = (real)NUM2DBL(rb_sin);

  sisnan_(&__out__, &sin);

  rb___out__ = __out__ ? Qtrue : Qfalse;
  return rb___out__;
}

void
init_lapack_sisnan(VALUE mLapack){
  rb_define_module_function(mLapack, "sisnan", rb_sisnan, -1);
}
