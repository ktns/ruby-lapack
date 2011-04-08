#include "rb_lapack.h"

extern VOID xerbla_(char *srname, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_xerbla(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_srname;
  char *srname; 
  VALUE rblapack_info;
  integer info; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n   = NumRu::Lapack.xerbla( srname, info, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE XERBLA( SRNAME, INFO )\n\n*  Purpose\n*  =======\n*\n*  XERBLA  is an error handler for the LAPACK routines.\n*  It is called by an LAPACK routine if an input parameter has an\n*  invalid value.  A message is printed and execution stops.\n*\n*  Installers may consider modifying the STOP statement in order to\n*  call system-specific exception-handling facilities.\n*\n\n*  Arguments\n*  =========\n*\n*  SRNAME  (input) CHARACTER*(*)\n*          The name of the routine which called XERBLA.\n*\n*  INFO    (input) INTEGER\n*          The position of the invalid parameter in the parameter list\n*          of the calling routine.\n*\n\n* =====================================================================\n*\n*     .. Intrinsic Functions ..\n      INTRINSIC          LEN_TRIM\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n   = NumRu::Lapack.xerbla( srname, info, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_srname = argv[0];
  rblapack_info = argv[1];
  if (rb_options != Qnil) {
  }

  srname = StringValueCStr(rblapack_srname);
  info = NUM2INT(rblapack_info);

  xerbla_(srname, &info);

  return Qnil;
}

void
init_lapack_xerbla(VALUE mLapack){
  rb_define_module_function(mLapack, "xerbla", rblapack_xerbla, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
