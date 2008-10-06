#include "rb_lapack.h"

static VALUE
rb_dlartg(int argc, VALUE *argv, VALUE self){
  VALUE rb_f;
  doublereal f; 
  VALUE rb_g;
  doublereal g; 
  VALUE rb_cs;
  doublereal cs; 
  VALUE rb_sn;
  doublereal sn; 
  VALUE rb_r;
  doublereal r; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.dlartg( f, g)\n    or\n  NumRu::Lapack.dlartg  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARTG( F, G, CS, SN, R )\n\n*  Purpose\n*  =======\n*\n*  DLARTG generate a plane rotation so that\n*\n*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.\n*     [ -SN  CS  ]     [ G ]     [ 0 ]\n*\n*  This is a slower, more accurate version of the BLAS1 routine DROTG,\n*  with the following other differences:\n*     F and G are unchanged on return.\n*     If G=0, then CS=1 and SN=0.\n*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any\n*        floating point operations (saves work in DBDSQR when\n*        there are zeros on the diagonal).\n*\n*  If F exceeds G in magnitude, CS will be positive.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) DOUBLE PRECISION\n*          The first component of vector to be rotated.\n*\n*  G       (input) DOUBLE PRECISION\n*          The second component of vector to be rotated.\n*\n*  CS      (output) DOUBLE PRECISION\n*          The cosine of the rotation.\n*\n*  SN      (output) DOUBLE PRECISION\n*          The sine of the rotation.\n*\n*  R       (output) DOUBLE PRECISION\n*          The nonzero component of the rotated vector.\n*\n*  This version has a few statements commented out for thread safety\n*  (machine parameters are computed on each entry). 10 feb 03, SJH.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_f = argv[0];
  rb_g = argv[1];

  f = NUM2DBL(rb_f);
  g = NUM2DBL(rb_g);

  dlartg_(&f, &g, &cs, &sn, &r);

  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_float_new((double)sn);
  rb_r = rb_float_new((double)r);
  return rb_ary_new3(3, rb_cs, rb_sn, rb_r);
}

void
init_lapack_dlartg(VALUE mLapack){
  rb_define_module_function(mLapack, "dlartg", rb_dlartg, -1);
}
