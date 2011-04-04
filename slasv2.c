#include "rb_lapack.h"

extern VOID slasv2_(real *f, real *g, real *h, real *ssmin, real *ssmax, real *snr, real *csr, real *snl, real *csl);

static VALUE
rb_slasv2(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_snr;
  real snr; 
  VALUE rb_csr;
  real csr; 
  VALUE rb_snl;
  real snl; 
  VALUE rb_csl;
  real csl; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ssmin, ssmax, snr, csr, snl, csl = NumRu::Lapack.slasv2( f, g, h)\n    or\n  NumRu::Lapack.slasv2  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  slasv2_(&f, &g, &h, &ssmin, &ssmax, &snr, &csr, &snl, &csl);

  rb_ssmin = rb_float_new((double)ssmin);
  rb_ssmax = rb_float_new((double)ssmax);
  rb_snr = rb_float_new((double)snr);
  rb_csr = rb_float_new((double)csr);
  rb_snl = rb_float_new((double)snl);
  rb_csl = rb_float_new((double)csl);
  return rb_ary_new3(6, rb_ssmin, rb_ssmax, rb_snr, rb_csr, rb_snl, rb_csl);
}

void
init_lapack_slasv2(VALUE mLapack){
  rb_define_module_function(mLapack, "slasv2", rb_slasv2, -1);
}
