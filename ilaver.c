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
    printf("%s\n", "USAGE:\n  vers_major, vers_minor, vers_patch = NumRu::Lapack.ilaver( )\n    or\n  NumRu::Lapack.ilaver  # print help\n\n\nFORTRAN MANUAL\n\n");
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
