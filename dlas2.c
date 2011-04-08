#include "rb_lapack.h"

extern VOID dlas2_(doublereal *f, doublereal *g, doublereal *h, doublereal *ssmin, doublereal *ssmax);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlas2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_f;
  doublereal f; 
  VALUE rblapack_g;
  doublereal g; 
  VALUE rblapack_h;
  doublereal h; 
  VALUE rblapack_ssmin;
  doublereal ssmin; 
  VALUE rblapack_ssmax;
  doublereal ssmax; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ssmin, ssmax = NumRu::Lapack.dlas2( f, g, h, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX )\n\n*  Purpose\n*  =======\n*\n*  DLAS2  computes the singular values of the 2-by-2 matrix\n*     [  F   G  ]\n*     [  0   H  ].\n*  On return, SSMIN is the smaller singular value and SSMAX is the\n*  larger singular value.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) DOUBLE PRECISION\n*          The (1,1) element of the 2-by-2 matrix.\n*\n*  G       (input) DOUBLE PRECISION\n*          The (1,2) element of the 2-by-2 matrix.\n*\n*  H       (input) DOUBLE PRECISION\n*          The (2,2) element of the 2-by-2 matrix.\n*\n*  SSMIN   (output) DOUBLE PRECISION\n*          The smaller singular value.\n*\n*  SSMAX   (output) DOUBLE PRECISION\n*          The larger singular value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Barring over/underflow, all output quantities are correct to within\n*  a few units in the last place (ulps), even in the absence of a guard\n*  digit in addition/subtraction.\n*\n*  In IEEE arithmetic, the code works correctly if one matrix element is\n*  infinite.\n*\n*  Overflow will not occur unless the largest singular value itself\n*  overflows, or is within a few ulps of overflow. (On machines with\n*  partial overflow, like the Cray, overflow may occur if the largest\n*  singular value is within a factor of 2 of overflow.)\n*\n*  Underflow is harmless if underflow is gradual. Otherwise, results\n*  may correspond to a matrix modified by perturbations of size near\n*  the underflow threshold.\n*\n*  ====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ssmin, ssmax = NumRu::Lapack.dlas2( f, g, h, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_f = argv[0];
  rblapack_g = argv[1];
  rblapack_h = argv[2];
  if (rb_options != Qnil) {
  }

  f = NUM2DBL(rblapack_f);
  g = NUM2DBL(rblapack_g);
  h = NUM2DBL(rblapack_h);

  dlas2_(&f, &g, &h, &ssmin, &ssmax);

  rblapack_ssmin = rb_float_new((double)ssmin);
  rblapack_ssmax = rb_float_new((double)ssmax);
  return rb_ary_new3(2, rblapack_ssmin, rblapack_ssmax);
}

void
init_lapack_dlas2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlas2", rblapack_dlas2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
