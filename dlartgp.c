#include "rb_lapack.h"

extern VOID dlartgp_(doublereal *f, doublereal *g, doublereal *cs, doublereal *sn, doublereal *r);

static VALUE
rb_dlartgp(int argc, VALUE *argv, VALUE self){
  VALUE rb_f;
  doublereal f; 
  VALUE rb_g;
  doublereal g; 
  VALUE rb_cs;
  doublereal cs; 
  VALUE rb_sn;
  doublereal sn; 
  VALUE rb_r;
  doublereal r; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.dlartgp( f, g)\n    or\n  NumRu::Lapack.dlartgp  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_f = argv[0];
  rb_g = argv[1];

  f = NUM2DBL(rb_f);
  g = NUM2DBL(rb_g);

  dlartgp_(&f, &g, &cs, &sn, &r);

  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_float_new((double)sn);
  rb_r = rb_float_new((double)r);
  return rb_ary_new3(3, rb_cs, rb_sn, rb_r);
}

void
init_lapack_dlartgp(VALUE mLapack){
  rb_define_module_function(mLapack, "dlartgp", rb_dlartgp, -1);
}
