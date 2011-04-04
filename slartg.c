#include "rb_lapack.h"

extern VOID slartg_(real *f, real *g, real *cs, real *sn, real *r);

static VALUE
rb_slartg(int argc, VALUE *argv, VALUE self){
  VALUE rb_f;
  real f; 
  VALUE rb_g;
  real g; 
  VALUE rb_cs;
  real cs; 
  VALUE rb_sn;
  real sn; 
  VALUE rb_r;
  real r; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.slartg( f, g)\n    or\n  NumRu::Lapack.slartg  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_f = argv[0];
  rb_g = argv[1];

  f = (real)NUM2DBL(rb_f);
  g = (real)NUM2DBL(rb_g);

  slartg_(&f, &g, &cs, &sn, &r);

  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_float_new((double)sn);
  rb_r = rb_float_new((double)r);
  return rb_ary_new3(3, rb_cs, rb_sn, rb_r);
}

void
init_lapack_slartg(VALUE mLapack){
  rb_define_module_function(mLapack, "slartg", rb_slartg, -1);
}
