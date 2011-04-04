#include "rb_lapack.h"

extern VOID dlae2_(doublereal *a, doublereal *b, doublereal *c, doublereal *rt1, doublereal *rt2);

static VALUE
rb_dlae2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal a; 
  VALUE rb_b;
  doublereal b; 
  VALUE rb_c;
  doublereal c; 
  VALUE rb_rt1;
  doublereal rt1; 
  VALUE rb_rt2;
  doublereal rt2; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rt1, rt2 = NumRu::Lapack.dlae2( a, b, c)\n    or\n  NumRu::Lapack.dlae2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_c = argv[2];

  a = NUM2DBL(rb_a);
  b = NUM2DBL(rb_b);
  c = NUM2DBL(rb_c);

  dlae2_(&a, &b, &c, &rt1, &rt2);

  rb_rt1 = rb_float_new((double)rt1);
  rb_rt2 = rb_float_new((double)rt2);
  return rb_ary_new3(2, rb_rt1, rb_rt2);
}

void
init_lapack_dlae2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlae2", rb_dlae2, -1);
}
