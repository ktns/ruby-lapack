#include "rb_lapack.h"

extern VOID dlartgs_(doublereal *x, doublereal *y, doublereal *sigma, doublereal *cs, doublereal *sn);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlartgs(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  doublereal x; 
  VALUE rblapack_y;
  doublereal y; 
  VALUE rblapack_sigma;
  doublereal sigma; 
  VALUE rblapack_cs;
  doublereal cs; 
  VALUE rblapack_sn;
  doublereal sn; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn = NumRu::Lapack.dlartgs( x, y, sigma, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARTGS( X, Y, SIGMA, CS, SN )\n\n*  Purpose\n*  =======\n*\n*  DLARTGS generates a plane rotation designed to introduce a bulge in\n*  Golub-Reinsch-style implicit QR iteration for the bidiagonal SVD\n*  problem. X and Y are the top-row entries, and SIGMA is the shift.\n*  The computed CS and SN define a plane rotation satisfying\n*\n*     [  CS  SN  ]  .  [ X^2 - SIGMA ]  =  [ R ],\n*     [ -SN  CS  ]     [    X * Y    ]     [ 0 ]\n*\n*  with R nonnegative.  If X^2 - SIGMA and X * Y are 0, then the\n*  rotation is by PI/2.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) DOUBLE PRECISION\n*          The (1,1) entry of an upper bidiagonal matrix.\n*\n*  Y       (input) DOUBLE PRECISION\n*          The (1,2) entry of an upper bidiagonal matrix.\n*\n*  SIGMA   (input) DOUBLE PRECISION\n*          The shift.\n*\n*  CS      (output) DOUBLE PRECISION\n*          The cosine of the rotation.\n*\n*  SN      (output) DOUBLE PRECISION\n*          The sine of the rotation.\n*\n\n*  ===================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn = NumRu::Lapack.dlartgs( x, y, sigma, [:usage => usage, :help => help])\n");
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

  x = NUM2DBL(rblapack_x);
  y = NUM2DBL(rblapack_y);
  sigma = NUM2DBL(rblapack_sigma);

  dlartgs_(&x, &y, &sigma, &cs, &sn);

  rblapack_cs = rb_float_new((double)cs);
  rblapack_sn = rb_float_new((double)sn);
  return rb_ary_new3(2, rblapack_cs, rblapack_sn);
}

void
init_lapack_dlartgs(VALUE mLapack){
  rb_define_module_function(mLapack, "dlartgs", rblapack_dlartgs, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
