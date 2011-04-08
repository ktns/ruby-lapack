#include "rb_lapack.h"

extern VOID slartgs_(real *x, real *y, real *sigma, real *cs, real *sn);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slartgs(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  real x; 
  VALUE rblapack_y;
  real y; 
  VALUE rblapack_sigma;
  real sigma; 
  VALUE rblapack_cs;
  real cs; 
  VALUE rblapack_sn;
  real sn; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn = NumRu::Lapack.slartgs( x, y, sigma, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARTGS( X, Y, SIGMA, CS, SN )\n\n*  Purpose\n*  =======\n*\n*  SLARTGS generates a plane rotation designed to introduce a bulge in\n*  Golub-Reinsch-style implicit QR iteration for the bidiagonal SVD\n*  problem. X and Y are the top-row entries, and SIGMA is the shift.\n*  The computed CS and SN define a plane rotation satisfying\n*\n*     [  CS  SN  ]  .  [ X^2 - SIGMA ]  =  [ R ],\n*     [ -SN  CS  ]     [    X * Y    ]     [ 0 ]\n*\n*  with R nonnegative.  If X^2 - SIGMA and X * Y are 0, then the\n*  rotation is by PI/2.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) REAL\n*          The (1,1) entry of an upper bidiagonal matrix.\n*\n*  Y       (input) REAL\n*          The (1,2) entry of an upper bidiagonal matrix.\n*\n*  SIGMA   (input) REAL\n*          The shift.\n*\n*  CS      (output) REAL\n*          The cosine of the rotation.\n*\n*  SN      (output) REAL\n*          The sine of the rotation.\n*\n\n*  ===================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn = NumRu::Lapack.slartgs( x, y, sigma, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_x = argv[0];
  rblapack_y = argv[1];
  rblapack_sigma = argv[2];
  if (rb_options != Qnil) {
  }

  x = (real)NUM2DBL(rblapack_x);
  y = (real)NUM2DBL(rblapack_y);
  sigma = (real)NUM2DBL(rblapack_sigma);

  slartgs_(&x, &y, &sigma, &cs, &sn);

  rblapack_cs = rb_float_new((double)cs);
  rblapack_sn = rb_float_new((double)sn);
  return rb_ary_new3(2, rblapack_cs, rblapack_sn);
}

void
init_lapack_slartgs(VALUE mLapack){
  rb_define_module_function(mLapack, "slartgs", rblapack_slartgs, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
