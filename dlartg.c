#include "rb_lapack.h"

extern VOID dlartg_(doublereal *f, doublereal *g, doublereal *cs, doublereal *sn, doublereal *r);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlartg(int argc, VALUE *argv, VALUE self){
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
      printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.dlartg( f, g, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARTG( F, G, CS, SN, R )\n\n*  Purpose\n*  =======\n*\n*  DLARTG generate a plane rotation so that\n*\n*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.\n*     [ -SN  CS  ]     [ G ]     [ 0 ]\n*\n*  This is a slower, more accurate version of the BLAS1 routine DROTG,\n*  with the following other differences:\n*     F and G are unchanged on return.\n*     If G=0, then CS=1 and SN=0.\n*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any\n*        floating point operations (saves work in DBDSQR when\n*        there are zeros on the diagonal).\n*\n*  If F exceeds G in magnitude, CS will be positive.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) DOUBLE PRECISION\n*          The first component of vector to be rotated.\n*\n*  G       (input) DOUBLE PRECISION\n*          The second component of vector to be rotated.\n*\n*  CS      (output) DOUBLE PRECISION\n*          The cosine of the rotation.\n*\n*  SN      (output) DOUBLE PRECISION\n*          The sine of the rotation.\n*\n*  R       (output) DOUBLE PRECISION\n*          The nonzero component of the rotated vector.\n*\n*  This version has a few statements commented out for thread safety\n*  (machine parameters are computed on each entry). 10 feb 03, SJH.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.dlartg( f, g, [:usage => usage, :help => help])\n");
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

  dlartg_(&f, &g, &cs, &sn, &r);

  rblapack_cs = rb_float_new((double)cs);
  rblapack_sn = rb_float_new((double)sn);
  rblapack_r = rb_float_new((double)r);
  return rb_ary_new3(3, rblapack_cs, rblapack_sn, rblapack_r);
}

void
init_lapack_dlartg(VALUE mLapack){
  rb_define_module_function(mLapack, "dlartg", rblapack_dlartg, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
