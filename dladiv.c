#include "rb_lapack.h"

extern VOID dladiv_(doublereal *a, doublereal *b, doublereal *c, doublereal *d, doublereal *p, doublereal *q);

static VALUE
rb_dladiv(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal a; 
  VALUE rb_b;
  doublereal b; 
  VALUE rb_c;
  doublereal c; 
  VALUE rb_d;
  doublereal d; 
  VALUE rb_p;
  doublereal p; 
  VALUE rb_q;
  doublereal q; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  p, q = NumRu::Lapack.dladiv( a, b, c, d)\n    or\n  NumRu::Lapack.dladiv  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  dladiv_(&a, &b, &c, &d, &p, &q);

  rb_p = rb_float_new((double)p);
  rb_q = rb_float_new((double)q);
  return rb_ary_new3(2, rb_p, rb_q);
}

void
init_lapack_dladiv(VALUE mLapack){
  rb_define_module_function(mLapack, "dladiv", rb_dladiv, -1);
}
