#include "rb_lapack.h"

extern VOID dlasv2_(doublereal *f, doublereal *g, doublereal *h, doublereal *ssmin, doublereal *ssmax, doublereal *snr, doublereal *csr, doublereal *snl, doublereal *csl);

static VALUE
rb_dlasv2(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_snr;
  doublereal snr; 
  VALUE rb_csr;
  doublereal csr; 
  VALUE rb_snl;
  doublereal snl; 
  VALUE rb_csl;
  doublereal csl; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ssmin, ssmax, snr, csr, snl, csl = NumRu::Lapack.dlasv2( f, g, h)\n    or\n  NumRu::Lapack.dlasv2  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  dlasv2_(&f, &g, &h, &ssmin, &ssmax, &snr, &csr, &snl, &csl);

  rb_ssmin = rb_float_new((double)ssmin);
  rb_ssmax = rb_float_new((double)ssmax);
  rb_snr = rb_float_new((double)snr);
  rb_csr = rb_float_new((double)csr);
  rb_snl = rb_float_new((double)snl);
  rb_csl = rb_float_new((double)csl);
  return rb_ary_new3(6, rb_ssmin, rb_ssmax, rb_snr, rb_csr, rb_snl, rb_csl);
}

void
init_lapack_dlasv2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasv2", rb_dlasv2, -1);
}
