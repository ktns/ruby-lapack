#include "rb_lapack.h"

extern VOID clartg_(complex *f, complex *g, real *cs, complex *sn, complex *r);

static VALUE
rb_clartg(int argc, VALUE *argv, VALUE self){
  VALUE rb_f;
  complex f; 
  VALUE rb_g;
  complex g; 
  VALUE rb_cs;
  real cs; 
  VALUE rb_sn;
  complex sn; 
  VALUE rb_r;
  complex r; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.clartg( f, g)\n    or\n  NumRu::Lapack.clartg  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_f = argv[0];
  rb_g = argv[1];

  f.r = (real)NUM2DBL(rb_funcall(rb_f, rb_intern("real"), 0));
  f.i = (real)NUM2DBL(rb_funcall(rb_f, rb_intern("imag"), 0));
  g.r = (real)NUM2DBL(rb_funcall(rb_g, rb_intern("real"), 0));
  g.i = (real)NUM2DBL(rb_funcall(rb_g, rb_intern("imag"), 0));

  clartg_(&f, &g, &cs, &sn, &r);

  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(sn.r)), rb_float_new((double)(sn.i)));
  rb_r = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(r.r)), rb_float_new((double)(r.i)));
  return rb_ary_new3(3, rb_cs, rb_sn, rb_r);
}

void
init_lapack_clartg(VALUE mLapack){
  rb_define_module_function(mLapack, "clartg", rb_clartg, -1);
}
