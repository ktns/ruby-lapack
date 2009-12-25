#include "rb_lapack.h"

static VALUE
rb_lsamen(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_ca;
  char *ca; 
  VALUE rb_cb;
  char *cb; 
  VALUE rb___out__;
  logical __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.lsamen( n, ca, cb)\n    or\n  NumRu::Lapack.lsamen  # print help\n\n\nFORTRAN MANUAL\n      LOGICAL          FUNCTION LSAMEN( N, CA, CB )\n\n*  Purpose\n*  =======\n*\n*  LSAMEN  tests if the first N letters of CA are the same as the\n*  first N letters of CB, regardless of case.\n*  LSAMEN returns .TRUE. if CA and CB are equivalent except for case\n*  and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA )\n*  or LEN( CB ) is less than N.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of characters in CA and CB to be compared.\n*\n*  CA      (input) CHARACTER*(*)\n*  CB      (input) CHARACTER*(*)\n*          CA and CB specify two character strings of length at least N.\n*          Only the first N characters of each string will be accessed.\n*\n\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          LEN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_n = argv[0];
  rb_ca = argv[1];
  rb_cb = argv[2];

  n = NUM2INT(rb_n);
  ca = StringValueCStr(rb_ca);
  cb = StringValueCStr(rb_cb);

  __out__ = lsamen_(&n, ca, cb);

  rb___out__ = __out__ ? Qtrue : Qfalse;
  return rb___out__;
}

void
init_lapack_lsamen(VALUE mLapack){
  rb_define_module_function(mLapack, "lsamen", rb_lsamen, -1);
}
