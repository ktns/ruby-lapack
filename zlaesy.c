#include "rb_lapack.h"

extern VOID zlaesy_(doublecomplex *a, doublecomplex *b, doublecomplex *c, doublecomplex *rt1, doublecomplex *rt2, doublecomplex *evscal, doublecomplex *cs1, doublecomplex *sn1);

static VALUE
rb_zlaesy(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublecomplex a; 
  VALUE rb_b;
  doublecomplex b; 
  VALUE rb_c;
  doublecomplex c; 
  VALUE rb_rt1;
  doublecomplex rt1; 
  VALUE rb_rt2;
  doublecomplex rt2; 
  VALUE rb_evscal;
  doublecomplex evscal; 
  VALUE rb_cs1;
  doublecomplex cs1; 
  VALUE rb_sn1;
  doublecomplex sn1; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rt1, rt2, evscal, cs1, sn1 = NumRu::Lapack.zlaesy( a, b, c)\n    or\n  NumRu::Lapack.zlaesy  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  zlaesy_(&a, &b, &c, &rt1, &rt2, &evscal, &cs1, &sn1);

  rb_rt1 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(rt1.r)), rb_float_new((double)(rt1.i)));
  rb_rt2 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(rt2.r)), rb_float_new((double)(rt2.i)));
  rb_evscal = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(evscal.r)), rb_float_new((double)(evscal.i)));
  rb_cs1 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(cs1.r)), rb_float_new((double)(cs1.i)));
  rb_sn1 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(sn1.r)), rb_float_new((double)(sn1.i)));
  return rb_ary_new3(5, rb_rt1, rb_rt2, rb_evscal, rb_cs1, rb_sn1);
}

void
init_lapack_zlaesy(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaesy", rb_zlaesy, -1);
}
