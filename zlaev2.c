#include "rb_lapack.h"

extern VOID zlaev2_(doublecomplex *a, doublecomplex *b, doublecomplex *c, doublereal *rt1, doublereal *rt2, doublereal *cs1, doublecomplex *sn1);

static VALUE
rb_zlaev2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublecomplex a; 
  VALUE rb_b;
  doublecomplex b; 
  VALUE rb_c;
  doublecomplex c; 
  VALUE rb_rt1;
  doublereal rt1; 
  VALUE rb_rt2;
  doublereal rt2; 
  VALUE rb_cs1;
  doublereal cs1; 
  VALUE rb_sn1;
  doublecomplex sn1; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rt1, rt2, cs1, sn1 = NumRu::Lapack.zlaev2( a, b, c)\n    or\n  NumRu::Lapack.zlaev2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_c = argv[2];

  a.r = NUM2DBL(rb_funcall(rb_a, rb_intern("real"), 0));
  a.i = NUM2DBL(rb_funcall(rb_a, rb_intern("imag"), 0));
  b.r = NUM2DBL(rb_funcall(rb_b, rb_intern("real"), 0));
  b.i = NUM2DBL(rb_funcall(rb_b, rb_intern("imag"), 0));
  c.r = NUM2DBL(rb_funcall(rb_c, rb_intern("real"), 0));
  c.i = NUM2DBL(rb_funcall(rb_c, rb_intern("imag"), 0));

  zlaev2_(&a, &b, &c, &rt1, &rt2, &cs1, &sn1);

  rb_rt1 = rb_float_new((double)rt1);
  rb_rt2 = rb_float_new((double)rt2);
  rb_cs1 = rb_float_new((double)cs1);
  rb_sn1 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(sn1.r)), rb_float_new((double)(sn1.i)));
  return rb_ary_new3(4, rb_rt1, rb_rt2, rb_cs1, rb_sn1);
}

void
init_lapack_zlaev2(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaev2", rb_zlaev2, -1);
}
