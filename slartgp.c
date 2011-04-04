#include "rb_lapack.h"

extern VOID slartgp_(real *f, real *g, real *cs, real *sn, real *r);

static VALUE
rb_slartgp(int argc, VALUE *argv, VALUE self){
  VALUE rb_f;
  real f; 
  VALUE rb_g;
  real g; 
  VALUE rb_cs;
  real cs; 
  VALUE rb_sn;
  real sn; 
  VALUE rb_r;
  real r; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.slartgp( f, g)\n    or\n  NumRu::Lapack.slartgp  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARTGP( F, G, CS, SN, R )\n\n*  Purpose\n*  =======\n*\n*  SLARTGP generates a plane rotation so that\n*\n*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.\n*     [ -SN  CS  ]     [ G ]     [ 0 ]\n*\n*  This is a slower, more accurate version of the Level 1 BLAS routine SROTG,\n*  with the following other differences:\n*     F and G are unchanged on return.\n*     If G=0, then CS=(+/-)1 and SN=0.\n*     If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1.\n*\n*  The sign is chosen so that R >= 0.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) REAL\n*          The first component of vector to be rotated.\n*\n*  G       (input) REAL\n*          The second component of vector to be rotated.\n*\n*  CS      (output) REAL\n*          The cosine of the rotation.\n*\n*  SN      (output) REAL\n*          The sine of the rotation.\n*\n*  R       (output) REAL\n*          The nonzero component of the rotated vector.\n*\n*  This version has a few statements commented out for thread safety\n*  (machine parameters are computed on each entry). 10 feb 03, SJH.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_f = argv[0];
  rb_g = argv[1];

  f = (real)NUM2DBL(rb_f);
  g = (real)NUM2DBL(rb_g);

  slartgp_(&f, &g, &cs, &sn, &r);

  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_float_new((double)sn);
  rb_r = rb_float_new((double)r);
  return rb_ary_new3(3, rb_cs, rb_sn, rb_r);
}

void
init_lapack_slartgp(VALUE mLapack){
  rb_define_module_function(mLapack, "slartgp", rb_slartgp, -1);
}
