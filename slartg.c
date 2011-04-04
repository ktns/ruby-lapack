#include "rb_lapack.h"

extern VOID slartg_(real *f, real *g, real *cs, real *sn, real *r);

static VALUE
rb_slartg(int argc, VALUE *argv, VALUE self){
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
    printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.slartg( f, g)\n    or\n  NumRu::Lapack.slartg  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARTG( F, G, CS, SN, R )\n\n*  Purpose\n*  =======\n*\n*  SLARTG generate a plane rotation so that\n*\n*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.\n*     [ -SN  CS  ]     [ G ]     [ 0 ]\n*\n*  This is a slower, more accurate version of the BLAS1 routine SROTG,\n*  with the following other differences:\n*     F and G are unchanged on return.\n*     If G=0, then CS=1 and SN=0.\n*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any\n*        floating point operations (saves work in SBDSQR when\n*        there are zeros on the diagonal).\n*\n*  If F exceeds G in magnitude, CS will be positive.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) REAL\n*          The first component of vector to be rotated.\n*\n*  G       (input) REAL\n*          The second component of vector to be rotated.\n*\n*  CS      (output) REAL\n*          The cosine of the rotation.\n*\n*  SN      (output) REAL\n*          The sine of the rotation.\n*\n*  R       (output) REAL\n*          The nonzero component of the rotated vector.\n*\n*  This version has a few statements commented out for thread safety\n*  (machine parameters are computed on each entry). 10 feb 03, SJH.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_f = argv[0];
  rb_g = argv[1];

  f = (real)NUM2DBL(rb_f);
  g = (real)NUM2DBL(rb_g);

  slartg_(&f, &g, &cs, &sn, &r);

  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_float_new((double)sn);
  rb_r = rb_float_new((double)r);
  return rb_ary_new3(3, rb_cs, rb_sn, rb_r);
}

void
init_lapack_slartg(VALUE mLapack){
  rb_define_module_function(mLapack, "slartg", rb_slartg, -1);
}
