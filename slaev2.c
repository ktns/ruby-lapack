#include "rb_lapack.h"

extern VOID slaev2_(real *a, real *b, real *c, real *rt1, real *rt2, real *cs1, real *sn1);

static VALUE
rb_slaev2(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_cs1;
  real cs1; 
  VALUE rb_sn1;
  real sn1; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rt1, rt2, cs1, sn1 = NumRu::Lapack.slaev2( a, b, c)\n    or\n  NumRu::Lapack.slaev2  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  slaev2_(&a, &b, &c, &rt1, &rt2, &cs1, &sn1);

  rb_rt1 = rb_float_new((double)rt1);
  rb_rt2 = rb_float_new((double)rt2);
  rb_cs1 = rb_float_new((double)cs1);
  rb_sn1 = rb_float_new((double)sn1);
  return rb_ary_new3(4, rb_rt1, rb_rt2, rb_cs1, rb_sn1);
}

void
init_lapack_slaev2(VALUE mLapack){
  rb_define_module_function(mLapack, "slaev2", rb_slaev2, -1);
}
