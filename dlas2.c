#include "rb_lapack.h"

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
    printf("%s\n", "USAGE:\n  ssmin, ssmax = NumRu::Lapack.dlas2( f, g, h)\n    or\n  NumRu::Lapack.dlas2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX )\n\n*  Purpose\n*  =======\n*\n*  DLAS2  computes the singular values of the 2-by-2 matrix\n*     [  F   G  ]\n*     [  0   H  ].\n*  On return, SSMIN is the smaller singular value and SSMAX is the\n*  larger singular value.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) DOUBLE PRECISION\n*          The (1,1) element of the 2-by-2 matrix.\n*\n*  G       (input) DOUBLE PRECISION\n*          The (1,2) element of the 2-by-2 matrix.\n*\n*  H       (input) DOUBLE PRECISION\n*          The (2,2) element of the 2-by-2 matrix.\n*\n*  SSMIN   (output) DOUBLE PRECISION\n*          The smaller singular value.\n*\n*  SSMAX   (output) DOUBLE PRECISION\n*          The larger singular value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Barring over/underflow, all output quantities are correct to within\n*  a few units in the last place (ulps), even in the absence of a guard\n*  digit in addition/subtraction.\n*\n*  In IEEE arithmetic, the code works correctly if one matrix element is\n*  infinite.\n*\n*  Overflow will not occur unless the largest singular value itself\n*  overflows, or is within a few ulps of overflow. (On machines with\n*  partial overflow, like the Cray, overflow may occur if the largest\n*  singular value is within a factor of 2 of overflow.)\n*\n*  Underflow is harmless if underflow is gradual. Otherwise, results\n*  may correspond to a matrix modified by perturbations of size near\n*  the underflow threshold.\n*\n*  ====================================================================\n*\n\n");
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
