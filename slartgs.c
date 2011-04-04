#include "rb_lapack.h"

extern VOID slartgs_(real *x, real *y, real *sigma, real *cs, real *sn);

static VALUE
rb_slartgs(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  real x; 
  VALUE rb_y;
  real y; 
  VALUE rb_sigma;
  real sigma; 
  VALUE rb_cs;
  real cs; 
  VALUE rb_sn;
  real sn; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  cs, sn = NumRu::Lapack.slartgs( x, y, sigma)\n    or\n  NumRu::Lapack.slartgs  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_x = argv[0];
  rb_y = argv[1];
  rb_sigma = argv[2];

  x = (real)NUM2DBL(rb_x);
  y = (real)NUM2DBL(rb_y);
  sigma = (real)NUM2DBL(rb_sigma);

  slartgs_(&x, &y, &sigma, &cs, &sn);

  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_float_new((double)sn);
  return rb_ary_new3(2, rb_cs, rb_sn);
}

void
init_lapack_slartgs(VALUE mLapack){
  rb_define_module_function(mLapack, "slartgs", rb_slartgs, -1);
}
