#include "rb_lapack.h"

extern VOID dlaev2_(doublereal *a, doublereal *b, doublereal *c, doublereal *rt1, doublereal *rt2, doublereal *cs1, doublereal *sn1);

static VALUE
rb_dlaev2(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_cs1;
  doublereal cs1; 
  VALUE rb_sn1;
  doublereal sn1; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rt1, rt2, cs1, sn1 = NumRu::Lapack.dlaev2( a, b, c)\n    or\n  NumRu::Lapack.dlaev2  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  dlaev2_(&a, &b, &c, &rt1, &rt2, &cs1, &sn1);

  rb_rt1 = rb_float_new((double)rt1);
  rb_rt2 = rb_float_new((double)rt2);
  rb_cs1 = rb_float_new((double)cs1);
  rb_sn1 = rb_float_new((double)sn1);
  return rb_ary_new3(4, rb_rt1, rb_rt2, rb_cs1, rb_sn1);
}

void
init_lapack_dlaev2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaev2", rb_dlaev2, -1);
}
