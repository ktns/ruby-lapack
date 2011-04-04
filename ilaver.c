#include "rb_lapack.h"

extern VOID ilaver_(integer *vers_major, integer *vers_minor, integer *vers_patch);

static VALUE
rb_ilaver(int argc, VALUE *argv, VALUE self){
  VALUE rb_vers_major;
  integer vers_major; 
  VALUE rb_vers_minor;
  integer vers_minor; 
  VALUE rb_vers_patch;
  integer vers_patch; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  vers_major, vers_minor, vers_patch = NumRu::Lapack.ilaver( )\n    or\n  NumRu::Lapack.ilaver  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )\n\n*  Purpose\n*  =======\n*\n*  This subroutine return the Lapack version.\n*\n\n*  Arguments\n*  =========\n*  VERS_MAJOR   (output) INTEGER\n*      return the lapack major version\n*  VERS_MINOR   (output) INTEGER\n*      return the lapack minor version from the major version\n*  VERS_PATCH   (output) INTEGER\n*      return the lapack patch version from the minor version\n\n*  =====================================================================\n*\n      INTEGER VERS_MAJOR, VERS_MINOR, VERS_PATCH\n*  =====================================================================\n      VERS_MAJOR = 3\n      VERS_MINOR = 3\n      VERS_PATCH = 0\n*  =====================================================================\n*\n      RETURN\n      END\n\n");
    return Qnil;
  }
  if (argc != 0)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 0)", argc);


  ilaver_(&vers_major, &vers_minor, &vers_patch);

  rb_vers_major = INT2NUM(vers_major);
  rb_vers_minor = INT2NUM(vers_minor);
  rb_vers_patch = INT2NUM(vers_patch);
  return rb_ary_new3(3, rb_vers_major, rb_vers_minor, rb_vers_patch);
}

void
init_lapack_ilaver(VALUE mLapack){
  rb_define_module_function(mLapack, "ilaver", rb_ilaver, -1);
}
