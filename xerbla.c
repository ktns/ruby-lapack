#include "rb_lapack.h"

extern VOID xerbla_(char *srname, integer *info);

static VALUE
rb_xerbla(int argc, VALUE *argv, VALUE self){
  VALUE rb_srname;
  char *srname; 
  VALUE rb_info;
  integer info; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n   = NumRu::Lapack.xerbla( srname, info)\n    or\n  NumRu::Lapack.xerbla  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE XERBLA( SRNAME, INFO )\n\n*  Purpose\n*  =======\n*\n*  XERBLA  is an error handler for the LAPACK routines.\n*  It is called by an LAPACK routine if an input parameter has an\n*  invalid value.  A message is printed and execution stops.\n*\n*  Installers may consider modifying the STOP statement in order to\n*  call system-specific exception-handling facilities.\n*\n\n*  Arguments\n*  =========\n*\n*  SRNAME  (input) CHARACTER*(*)\n*          The name of the routine which called XERBLA.\n*\n*  INFO    (input) INTEGER\n*          The position of the invalid parameter in the parameter list\n*          of the calling routine.\n*\n\n* =====================================================================\n*\n*     .. Intrinsic Functions ..\n      INTRINSIC          LEN_TRIM\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_srname = argv[0];
  rb_info = argv[1];

  srname = StringValueCStr(rb_srname);
  info = NUM2INT(rb_info);

  xerbla_(srname, &info);

  return Qnil;
}

void
init_lapack_xerbla(VALUE mLapack){
  rb_define_module_function(mLapack, "xerbla", rb_xerbla, -1);
}
