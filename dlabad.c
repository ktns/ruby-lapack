#include "rb_lapack.h"

extern VOID dlabad_(doublereal *small, doublereal *large);

static VALUE
rb_dlabad(int argc, VALUE *argv, VALUE self){
  VALUE rb_small;
  doublereal small; 
  VALUE rb_large;
  doublereal large; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  small, large = NumRu::Lapack.dlabad( small, large)\n    or\n  NumRu::Lapack.dlabad  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_small = argv[0];
  rb_large = argv[1];

  large = NUM2DBL(rb_large);
  small = NUM2DBL(rb_small);

  dlabad_(&small, &large);

  rb_small = rb_float_new((double)small);
  rb_large = rb_float_new((double)large);
  return rb_ary_new3(2, rb_small, rb_large);
}

void
init_lapack_dlabad(VALUE mLapack){
  rb_define_module_function(mLapack, "dlabad", rb_dlabad, -1);
}
