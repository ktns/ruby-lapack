#include "rb_lapack.h"

extern VOID xerbla_array_(char *srname_array, integer *srname_len, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_xerbla_array(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_srname_array;
  char *srname_array; 
  VALUE rblapack_info;
  integer info; 

  integer srname_len;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n   = NumRu::Lapack.xerbla_array( srname_array, info, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE XERBLA_ARRAY( SRNAME_ARRAY, SRNAME_LEN, INFO)\n\n*  Purpose\n*  =======\n*\n*  XERBLA_ARRAY assists other languages in calling XERBLA, the LAPACK\n*  and BLAS error handler.  Rather than taking a Fortran string argument\n*  as the function's name, XERBLA_ARRAY takes an array of single\n*  characters along with the array's length.  XERBLA_ARRAY then copies\n*  up to 32 characters of that array into a Fortran string and passes\n*  that to XERBLA.  If called with a non-positive SRNAME_LEN,\n*  XERBLA_ARRAY will call XERBLA with a string of all blank characters.\n*\n*  Say some macro or other device makes XERBLA_ARRAY available to C99\n*  by a name lapack_xerbla and with a common Fortran calling convention.\n*  Then a C99 program could invoke XERBLA via:\n*     {\n*       int flen = strlen(__func__);\n*       lapack_xerbla(__func__, &flen, &info);\n*     }\n*\n*  Providing XERBLA_ARRAY is not necessary for intercepting LAPACK\n*  errors.  XERBLA_ARRAY calls XERBLA.\n*\n\n*  Arguments\n*  =========\n*\n*  SRNAME_ARRAY (input) CHARACTER(1) array, dimension (SRNAME_LEN)\n*          The name of the routine which called XERBLA_ARRAY.\n*\n*  SRNAME_LEN (input) INTEGER\n*          The length of the name in SRNAME_ARRAY.\n*\n*  INFO    (input) INTEGER\n*          The position of the invalid parameter in the parameter list\n*          of the calling routine.\n*\n\n* =====================================================================\n*\n*     ..\n*     .. Local Scalars ..\n      INTEGER I\n*     ..\n*     .. Local Arrays ..\n      CHARACTER*32 SRNAME\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC MIN, LEN\n*     ..\n*     .. External Functions ..\n      EXTERNAL XERBLA\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n   = NumRu::Lapack.xerbla_array( srname_array, info, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_srname_array = argv[0];
  rblapack_info = argv[1];
  if (rb_options != Qnil) {
  }

  info = NUM2INT(rblapack_info);
  srname_array = StringValueCStr(rblapack_srname_array);

  xerbla_array_(srname_array, &srname_len, &info);

  return Qnil;
}

void
init_lapack_xerbla_array(VALUE mLapack){
  rb_define_module_function(mLapack, "xerbla_array", rblapack_xerbla_array, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
