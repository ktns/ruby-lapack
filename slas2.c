#include "rb_lapack.h"

extern VOID slas2_(real *f, real *g, real *h, real *ssmin, real *ssmax);

static VALUE
rb_slas2(int argc, VALUE *argv, VALUE self){
  VALUE rb_f;
  real f; 
  VALUE rb_g;
  real g; 
  VALUE rb_h;
  real h; 
  VALUE rb_ssmin;
  real ssmin; 
  VALUE rb_ssmax;
  real ssmax; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ssmin, ssmax = NumRu::Lapack.slas2( f, g, h)\n    or\n  NumRu::Lapack.slas2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_f = argv[0];
  rb_g = argv[1];
  rb_h = argv[2];

  f = (real)NUM2DBL(rb_f);
  g = (real)NUM2DBL(rb_g);
  h = (real)NUM2DBL(rb_h);

  slas2_(&f, &g, &h, &ssmin, &ssmax);

  rb_ssmin = rb_float_new((double)ssmin);
  rb_ssmax = rb_float_new((double)ssmax);
  return rb_ary_new3(2, rb_ssmin, rb_ssmax);
}

void
init_lapack_slas2(VALUE mLapack){
  rb_define_module_function(mLapack, "slas2", rb_slas2, -1);
}
