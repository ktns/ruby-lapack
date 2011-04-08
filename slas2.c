#include "rb_lapack.h"

extern VOID slas2_(real *f, real *g, real *h, real *ssmin, real *ssmax);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slas2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_f;
  real f; 
  VALUE rblapack_g;
  real g; 
  VALUE rblapack_h;
  real h; 
  VALUE rblapack_ssmin;
  real ssmin; 
  VALUE rblapack_ssmax;
  real ssmax; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ssmin, ssmax = NumRu::Lapack.slas2( f, g, h, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAS2( F, G, H, SSMIN, SSMAX )\n\n*  Purpose\n*  =======\n*\n*  SLAS2  computes the singular values of the 2-by-2 matrix\n*     [  F   G  ]\n*     [  0   H  ].\n*  On return, SSMIN is the smaller singular value and SSMAX is the\n*  larger singular value.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) REAL\n*          The (1,1) element of the 2-by-2 matrix.\n*\n*  G       (input) REAL\n*          The (1,2) element of the 2-by-2 matrix.\n*\n*  H       (input) REAL\n*          The (2,2) element of the 2-by-2 matrix.\n*\n*  SSMIN   (output) REAL\n*          The smaller singular value.\n*\n*  SSMAX   (output) REAL\n*          The larger singular value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Barring over/underflow, all output quantities are correct to within\n*  a few units in the last place (ulps), even in the absence of a guard\n*  digit in addition/subtraction.\n*\n*  In IEEE arithmetic, the code works correctly if one matrix element is\n*  infinite.\n*\n*  Overflow will not occur unless the largest singular value itself\n*  overflows, or is within a few ulps of overflow. (On machines with\n*  partial overflow, like the Cray, overflow may occur if the largest\n*  singular value is within a factor of 2 of overflow.)\n*\n*  Underflow is harmless if underflow is gradual. Otherwise, results\n*  may correspond to a matrix modified by perturbations of size near\n*  the underflow threshold.\n*\n*  ====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ssmin, ssmax = NumRu::Lapack.slas2( f, g, h, [:usage => usage, :help => help])\n");
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

  f = (real)NUM2DBL(rblapack_f);
  g = (real)NUM2DBL(rblapack_g);
  h = (real)NUM2DBL(rblapack_h);

  slas2_(&f, &g, &h, &ssmin, &ssmax);

  rblapack_ssmin = rb_float_new((double)ssmin);
  rblapack_ssmax = rb_float_new((double)ssmax);
  return rb_ary_new3(2, rblapack_ssmin, rblapack_ssmax);
}

void
init_lapack_slas2(VALUE mLapack){
  rb_define_module_function(mLapack, "slas2", rblapack_slas2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
