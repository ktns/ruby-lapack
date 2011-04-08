#include "rb_lapack.h"

extern VOID slartg_(real *f, real *g, real *cs, real *sn, real *r);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slartg(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_f;
  real f; 
  VALUE rblapack_g;
  real g; 
  VALUE rblapack_cs;
  real cs; 
  VALUE rblapack_sn;
  real sn; 
  VALUE rblapack_r;
  real r; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.slartg( f, g, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARTG( F, G, CS, SN, R )\n\n*  Purpose\n*  =======\n*\n*  SLARTG generate a plane rotation so that\n*\n*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.\n*     [ -SN  CS  ]     [ G ]     [ 0 ]\n*\n*  This is a slower, more accurate version of the BLAS1 routine SROTG,\n*  with the following other differences:\n*     F and G are unchanged on return.\n*     If G=0, then CS=1 and SN=0.\n*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any\n*        floating point operations (saves work in SBDSQR when\n*        there are zeros on the diagonal).\n*\n*  If F exceeds G in magnitude, CS will be positive.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) REAL\n*          The first component of vector to be rotated.\n*\n*  G       (input) REAL\n*          The second component of vector to be rotated.\n*\n*  CS      (output) REAL\n*          The cosine of the rotation.\n*\n*  SN      (output) REAL\n*          The sine of the rotation.\n*\n*  R       (output) REAL\n*          The nonzero component of the rotated vector.\n*\n*  This version has a few statements commented out for thread safety\n*  (machine parameters are computed on each entry). 10 feb 03, SJH.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.slartg( f, g, [:usage => usage, :help => help])\n");
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

  f = (real)NUM2DBL(rblapack_f);
  g = (real)NUM2DBL(rblapack_g);

  slartg_(&f, &g, &cs, &sn, &r);

  rblapack_cs = rb_float_new((double)cs);
  rblapack_sn = rb_float_new((double)sn);
  rblapack_r = rb_float_new((double)r);
  return rb_ary_new3(3, rblapack_cs, rblapack_sn, rblapack_r);
}

void
init_lapack_slartg(VALUE mLapack){
  rb_define_module_function(mLapack, "slartg", rblapack_slartg, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
