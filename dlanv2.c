#include "rb_lapack.h"

extern VOID dlanv2_(doublereal *a, doublereal *b, doublereal *c, doublereal *d, doublereal *rt1r, doublereal *rt1i, doublereal *rt2r, doublereal *rt2i, doublereal *cs, doublereal *sn);

static VALUE
rb_dlanv2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal a; 
  VALUE rb_b;
  doublereal b; 
  VALUE rb_c;
  doublereal c; 
  VALUE rb_d;
  doublereal d; 
  VALUE rb_rt1r;
  doublereal rt1r; 
  VALUE rb_rt1i;
  doublereal rt1i; 
  VALUE rb_rt2r;
  doublereal rt2r; 
  VALUE rb_rt2i;
  doublereal rt2i; 
  VALUE rb_cs;
  doublereal cs; 
  VALUE rb_sn;
  doublereal sn; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rt1r, rt1i, rt2r, rt2i, cs, sn, a, b, c, d = NumRu::Lapack.dlanv2( a, b, c, d)\n    or\n  NumRu::Lapack.dlanv2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_c = argv[2];
  rb_d = argv[3];

  a = NUM2DBL(rb_a);
  b = NUM2DBL(rb_b);
  c = NUM2DBL(rb_c);
  d = NUM2DBL(rb_d);

  dlanv2_(&a, &b, &c, &d, &rt1r, &rt1i, &rt2r, &rt2i, &cs, &sn);

  rb_rt1r = rb_float_new((double)rt1r);
  rb_rt1i = rb_float_new((double)rt1i);
  rb_rt2r = rb_float_new((double)rt2r);
  rb_rt2i = rb_float_new((double)rt2i);
  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_float_new((double)sn);
  rb_a = rb_float_new((double)a);
  rb_b = rb_float_new((double)b);
  rb_c = rb_float_new((double)c);
  rb_d = rb_float_new((double)d);
  return rb_ary_new3(10, rb_rt1r, rb_rt1i, rb_rt2r, rb_rt2i, rb_cs, rb_sn, rb_a, rb_b, rb_c, rb_d);
}

void
init_lapack_dlanv2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlanv2", rb_dlanv2, -1);
}
