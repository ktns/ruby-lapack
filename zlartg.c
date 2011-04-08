#include "rb_lapack.h"

extern VOID zlartg_(doublecomplex *f, doublecomplex *g, doublereal *cs, doublecomplex *sn, doublecomplex *r);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlartg(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_f;
  doublecomplex f; 
  VALUE rblapack_g;
  doublecomplex g; 
  VALUE rblapack_cs;
  doublereal cs; 
  VALUE rblapack_sn;
  doublecomplex sn; 
  VALUE rblapack_r;
  doublecomplex r; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.zlartg( f, g, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLARTG( F, G, CS, SN, R )\n\n*  Purpose\n*  =======\n*\n*  ZLARTG generates a plane rotation so that\n*\n*     [  CS  SN  ]     [ F ]     [ R ]\n*     [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.\n*     [ -SN  CS  ]     [ G ]     [ 0 ]\n*\n*  This is a faster version of the BLAS1 routine ZROTG, except for\n*  the following differences:\n*     F and G are unchanged on return.\n*     If G=0, then CS=1 and SN=0.\n*     If F=0, then CS=0 and SN is chosen so that R is real.\n*\n\n*  Arguments\n*  =========\n*\n*  F       (input) COMPLEX*16\n*          The first component of vector to be rotated.\n*\n*  G       (input) COMPLEX*16\n*          The second component of vector to be rotated.\n*\n*  CS      (output) DOUBLE PRECISION\n*          The cosine of the rotation.\n*\n*  SN      (output) COMPLEX*16\n*          The sine of the rotation.\n*\n*  R       (output) COMPLEX*16\n*          The nonzero component of the rotated vector.\n*\n\n*  Further Details\n*  ======= =======\n*\n*  3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel\n*\n*  This version has a few statements commented out for thread safety\n*  (machine parameters are computed on each entry). 10 feb 03, SJH.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  cs, sn, r = NumRu::Lapack.zlartg( f, g, [:usage => usage, :help => help])\n");
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

  f.r = NUM2DBL(rb_funcall(rblapack_f, rb_intern("real"), 0));
  f.i = NUM2DBL(rb_funcall(rblapack_f, rb_intern("imag"), 0));
  g.r = NUM2DBL(rb_funcall(rblapack_g, rb_intern("real"), 0));
  g.i = NUM2DBL(rb_funcall(rblapack_g, rb_intern("imag"), 0));

  zlartg_(&f, &g, &cs, &sn, &r);

  rblapack_cs = rb_float_new((double)cs);
  rblapack_sn = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(sn.r)), rb_float_new((double)(sn.i)));
  rblapack_r = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(r.r)), rb_float_new((double)(r.i)));
  return rb_ary_new3(3, rblapack_cs, rblapack_sn, rblapack_r);
}

void
init_lapack_zlartg(VALUE mLapack){
  rb_define_module_function(mLapack, "zlartg", rblapack_zlartg, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
