#include "rb_lapack.h"

static VALUE
rb_clartg(int argc, VALUE *argv, VALUE self){
  VALUE rb_f;
  complex f; 
  VALUE rb_g;
  complex g; 
  VALUE rb_cs;
  real cs; 
  VALUE rb_sn;
  complex sn; 
  VALUE rb_r;
  complex r; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.clartg( f, g)\n    or\n  NumRu::Lapack.clartg  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLARTG( F, G, CS, SN, R )\n\n*  Purpose\n*  =======\n*\n*  CLARTG generates a plane rotation so that\n*\n*     [  CS  SN  ]     [ F ]     [ R ]\n*     [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.\n*     [ -SN  CS  ]     [ G ]     [ 0 ]\n*\n*  This is a faster version of the BLAS1 routine CROTG, except for\n*  the following differences:\n*     F and G are unchanged on return.\n*     If G=0, then CS=1 and SN=0.\n*     If F=0, then CS=0 and SN is chosen so that R is real.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) COMPLEX\n*          The first component of vector to be rotated.\n*\n*  G       (input) COMPLEX\n*          The second component of vector to be rotated.\n*\n*  CS      (output) REAL\n*          The cosine of the rotation.\n*\n*  SN      (output) COMPLEX\n*          The sine of the rotation.\n*\n*  R       (output) COMPLEX\n*          The nonzero component of the rotated vector.\n*\n\n*  Further Details\n*  ======= =======\n*\n*  3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel\n*\n*  This version has a few statements commented out for thread safety\n*  (machine parameters are computed on each entry). 10 feb 03, SJH.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_f = argv[0];
  rb_g = argv[1];

  f.r = (real)NUM2DBL(rb_funcall(rb_f, rb_intern("real"), 0));
  f.i = (real)NUM2DBL(rb_funcall(rb_f, rb_intern("imag"), 0));
  g.r = (real)NUM2DBL(rb_funcall(rb_g, rb_intern("real"), 0));
  g.i = (real)NUM2DBL(rb_funcall(rb_g, rb_intern("imag"), 0));

  clartg_(&f, &g, &cs, &sn, &r);

  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(sn.r)), rb_float_new((double)(sn.i)));
  rb_r = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(r.r)), rb_float_new((double)(r.i)));
  return rb_ary_new3(3, rb_cs, rb_sn, rb_r);
}

void
init_lapack_clartg(VALUE mLapack){
  rb_define_module_function(mLapack, "clartg", rb_clartg, -1);
}
