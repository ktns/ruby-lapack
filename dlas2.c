#include "rb_lapack.h"

extern VOID dlas2_(doublereal *f, doublereal *g, doublereal *h, doublereal *ssmin, doublereal *ssmax);

static VALUE
rb_dlas2(int argc, VALUE *argv, VALUE self){
  VALUE rb_f;
  doublereal f; 
  VALUE rb_g;
  doublereal g; 
  VALUE rb_h;
  doublereal h; 
  VALUE rb_ssmin;
  doublereal ssmin; 
  VALUE rb_ssmax;
  doublereal ssmax; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ssmin, ssmax = NumRu::Lapack.dlas2( f, g, h)\n    or\n  NumRu::Lapack.dlas2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_f = argv[0];
  rb_g = argv[1];
  rb_h = argv[2];

  f = NUM2DBL(rb_f);
  g = NUM2DBL(rb_g);
  h = NUM2DBL(rb_h);

  dlas2_(&f, &g, &h, &ssmin, &ssmax);

  rb_ssmin = rb_float_new((double)ssmin);
  rb_ssmax = rb_float_new((double)ssmax);
  return rb_ary_new3(2, rb_ssmin, rb_ssmax);
}

void
init_lapack_dlas2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlas2", rb_dlas2, -1);
}
