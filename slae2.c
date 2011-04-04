#include "rb_lapack.h"

extern VOID slae2_(real *a, real *b, real *c, real *rt1, real *rt2);

static VALUE
rb_slae2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real a; 
  VALUE rb_b;
  real b; 
  VALUE rb_c;
  real c; 
  VALUE rb_rt1;
  real rt1; 
  VALUE rb_rt2;
  real rt2; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rt1, rt2 = NumRu::Lapack.slae2( a, b, c)\n    or\n  NumRu::Lapack.slae2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_c = argv[2];

  a = (real)NUM2DBL(rb_a);
  b = (real)NUM2DBL(rb_b);
  c = (real)NUM2DBL(rb_c);

  slae2_(&a, &b, &c, &rt1, &rt2);

  rb_rt1 = rb_float_new((double)rt1);
  rb_rt2 = rb_float_new((double)rt2);
  return rb_ary_new3(2, rb_rt1, rb_rt2);
}

void
init_lapack_slae2(VALUE mLapack){
  rb_define_module_function(mLapack, "slae2", rb_slae2, -1);
}
