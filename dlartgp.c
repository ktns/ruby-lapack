#include "rb_lapack.h"

extern VOID dlartgp_(doublereal *f, doublereal *g, doublereal *cs, doublereal *sn, doublereal *r);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlartgp(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_f;
  doublereal f; 
  VALUE rblapack_g;
  doublereal g; 
  VALUE rblapack_cs;
  doublereal cs; 
  VALUE rblapack_sn;
  doublereal sn; 
  VALUE rblapack_r;
  doublereal r; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.dlartgp( f, g, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARTGP( F, G, CS, SN, R )\n\n*  Purpose\n*  =======\n*\n*  DLARTGP generates a plane rotation so that\n*\n*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.\n*     [ -SN  CS  ]     [ G ]     [ 0 ]\n*\n*  This is a slower, more accurate version of the Level 1 BLAS routine DROTG,\n*  with the following other differences:\n*     F and G are unchanged on return.\n*     If G=0, then CS=(+/-)1 and SN=0.\n*     If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1.\n*\n*  The sign is chosen so that R >= 0.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) DOUBLE PRECISION\n*          The first component of vector to be rotated.\n*\n*  G       (input) DOUBLE PRECISION\n*          The second component of vector to be rotated.\n*\n*  CS      (output) DOUBLE PRECISION\n*          The cosine of the rotation.\n*\n*  SN      (output) DOUBLE PRECISION\n*          The sine of the rotation.\n*\n*  R       (output) DOUBLE PRECISION\n*          The nonzero component of the rotated vector.\n*\n*  This version has a few statements commented out for thread safety\n*  (machine parameters are computed on each entry). 10 feb 03, SJH.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.dlartgp( f, g, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_f = argv[0];
  rblapack_g = argv[1];
  if (rb_options != Qnil) {
  }

  f = NUM2DBL(rblapack_f);
  g = NUM2DBL(rblapack_g);

  dlartgp_(&f, &g, &cs, &sn, &r);

  rblapack_cs = rb_float_new((double)cs);
  rblapack_sn = rb_float_new((double)sn);
  rblapack_r = rb_float_new((double)r);
  return rb_ary_new3(3, rblapack_cs, rblapack_sn, rblapack_r);
}

void
init_lapack_dlartgp(VALUE mLapack){
  rb_define_module_function(mLapack, "dlartgp", rblapack_dlartgp, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
