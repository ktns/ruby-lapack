#include "rb_lapack.h"

extern logical lsamen_(integer *n, char *ca, char *cb);

static VALUE sHelp, sUsage;

static VALUE
rblapack_lsamen(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_ca;
  char *ca; 
  VALUE rblapack_cb;
  char *cb; 
  VALUE rblapack___out__;
  logical __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.lsamen( n, ca, cb, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      LOGICAL          FUNCTION LSAMEN( N, CA, CB )\n\n*  Purpose\n*  =======\n*\n*  LSAMEN  tests if the first N letters of CA are the same as the\n*  first N letters of CB, regardless of case.\n*  LSAMEN returns .TRUE. if CA and CB are equivalent except for case\n*  and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA )\n*  or LEN( CB ) is less than N.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of characters in CA and CB to be compared.\n*\n*  CA      (input) CHARACTER*(*)\n*  CB      (input) CHARACTER*(*)\n*          CA and CB specify two character strings of length at least N.\n*          Only the first N characters of each string will be accessed.\n*\n\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          LEN\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.lsamen( n, ca, cb, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_n = argv[0];
  rblapack_ca = argv[1];
  rblapack_cb = argv[2];
  if (rb_options != Qnil) {
  }

  n = NUM2INT(rblapack_n);
  ca = StringValueCStr(rblapack_ca);
  cb = StringValueCStr(rblapack_cb);

  __out__ = lsamen_(&n, ca, cb);

  rblapack___out__ = __out__ ? Qtrue : Qfalse;
  return rblapack___out__;
}

void
init_lapack_lsamen(VALUE mLapack){
  rb_define_module_function(mLapack, "lsamen", rblapack_lsamen, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
