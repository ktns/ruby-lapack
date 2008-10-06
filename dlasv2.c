#include "rb_lapack.h"

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
    printf("%s\n", "USAGE:\n  ssmin, ssmax, snr, csr, snl, csl = NumRu::Lapack.dlasv2( f, g, h)\n    or\n  NumRu::Lapack.dlasv2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )\n\n*  Purpose\n*  =======\n*\n*  DLASV2 computes the singular value decomposition of a 2-by-2\n*  triangular matrix\n*     [  F   G  ]\n*     [  0   H  ].\n*  On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the\n*  smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and\n*  right singular vectors for abs(SSMAX), giving the decomposition\n*\n*     [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]\n*     [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) DOUBLE PRECISION\n*          The (1,1) element of the 2-by-2 matrix.\n*\n*  G       (input) DOUBLE PRECISION\n*          The (1,2) element of the 2-by-2 matrix.\n*\n*  H       (input) DOUBLE PRECISION\n*          The (2,2) element of the 2-by-2 matrix.\n*\n*  SSMIN   (output) DOUBLE PRECISION\n*          abs(SSMIN) is the smaller singular value.\n*\n*  SSMAX   (output) DOUBLE PRECISION\n*          abs(SSMAX) is the larger singular value.\n*\n*  SNL     (output) DOUBLE PRECISION\n*  CSL     (output) DOUBLE PRECISION\n*          The vector (CSL, SNL) is a unit left singular vector for the\n*          singular value abs(SSMAX).\n*\n*  SNR     (output) DOUBLE PRECISION\n*  CSR     (output) DOUBLE PRECISION\n*          The vector (CSR, SNR) is a unit right singular vector for the\n*          singular value abs(SSMAX).\n*\n\n*  Further Details\n*  ===============\n*\n*  Any input parameter may be aliased with any output parameter.\n*\n*  Barring over/underflow and assuming a guard digit in subtraction, all\n*  output quantities are correct to within a few units in the last\n*  place (ulps).\n*\n*  In IEEE arithmetic, the code works correctly if one matrix element is\n*  infinite.\n*\n*  Overflow will not occur unless the largest singular value itself\n*  overflows or is within a few ulps of overflow. (On machines with\n*  partial overflow, like the Cray, overflow may occur if the largest\n*  singular value is within a factor of 2 of overflow.)\n*\n*  Underflow is harmless if underflow is gradual. Otherwise, results\n*  may correspond to a matrix modified by perturbations of size near\n*  the underflow threshold.\n*\n* =====================================================================\n*\n\n");
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
