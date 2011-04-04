#include "rb_lapack.h"

extern VOID dlartgs_(doublereal *x, doublereal *y, doublereal *sigma, doublereal *cs, doublereal *sn);

static VALUE
rb_dlartgs(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  doublereal x; 
  VALUE rb_y;
  doublereal y; 
  VALUE rb_sigma;
  doublereal sigma; 
  VALUE rb_cs;
  doublereal cs; 
  VALUE rb_sn;
  doublereal sn; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  cs, sn = NumRu::Lapack.dlartgs( x, y, sigma)\n    or\n  NumRu::Lapack.dlartgs  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_x = argv[0];
  rb_y = argv[1];
  rb_sigma = argv[2];

  x = NUM2DBL(rb_x);
  y = NUM2DBL(rb_y);
  sigma = NUM2DBL(rb_sigma);

  dlartgs_(&x, &y, &sigma, &cs, &sn);

  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_float_new((double)sn);
  return rb_ary_new3(2, rb_cs, rb_sn);
}

void
init_lapack_dlartgs(VALUE mLapack){
  rb_define_module_function(mLapack, "dlartgs", rb_dlartgs, -1);
}
