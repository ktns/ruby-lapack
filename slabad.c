#include "rb_lapack.h"

extern VOID slabad_(real *small, real *large);

static VALUE
rb_slabad(int argc, VALUE *argv, VALUE self){
  VALUE rb_small;
  real small; 
  VALUE rb_large;
  real large; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  small, large = NumRu::Lapack.slabad( small, large)\n    or\n  NumRu::Lapack.slabad  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_small = argv[0];
  rb_large = argv[1];

  large = (real)NUM2DBL(rb_large);
  small = (real)NUM2DBL(rb_small);

  slabad_(&small, &large);

  rb_small = rb_float_new((double)small);
  rb_large = rb_float_new((double)large);
  return rb_ary_new3(2, rb_small, rb_large);
}

void
init_lapack_slabad(VALUE mLapack){
  rb_define_module_function(mLapack, "slabad", rb_slabad, -1);
}
