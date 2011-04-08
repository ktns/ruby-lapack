#include "rb_lapack.h"

extern VOID slanv2_(real *a, real *b, real *c, real *d, real *rt1r, real *rt1i, real *rt2r, real *rt2i, real *cs, real *sn);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slanv2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  real a; 
  VALUE rblapack_b;
  real b; 
  VALUE rblapack_c;
  real c; 
  VALUE rblapack_d;
  real d; 
  VALUE rblapack_rt1r;
  real rt1r; 
  VALUE rblapack_rt1i;
  real rt1i; 
  VALUE rblapack_rt2r;
  real rt2r; 
  VALUE rblapack_rt2i;
  real rt2i; 
  VALUE rblapack_cs;
  real cs; 
  VALUE rblapack_sn;
  real sn; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rt1r, rt1i, rt2r, rt2i, cs, sn, a, b, c, d = NumRu::Lapack.slanv2( a, b, c, d, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )\n\n*  Purpose\n*  =======\n*\n*  SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric\n*  matrix in standard form:\n*\n*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]\n*       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]\n*\n*  where either\n*  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or\n*  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex\n*  conjugate eigenvalues.\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input/output) REAL            \n*  B       (input/output) REAL            \n*  C       (input/output) REAL            \n*  D       (input/output) REAL            \n*          On entry, the elements of the input matrix.\n*          On exit, they are overwritten by the elements of the\n*          standardised Schur form.\n*\n*  RT1R    (output) REAL \n*  RT1I    (output) REAL            \n*  RT2R    (output) REAL            \n*  RT2I    (output) REAL            \n*          The real and imaginary parts of the eigenvalues. If the\n*          eigenvalues are a complex conjugate pair, RT1I > 0.\n*\n*  CS      (output) REAL            \n*  SN      (output) REAL            \n*          Parameters of the rotation matrix.\n*\n\n*  Further Details\n*  ===============\n*\n*  Modified by V. Sima, Research Institute for Informatics, Bucharest,\n*  Romania, to reduce the risk of cancellation errors,\n*  when computing real eigenvalues, and to ensure, if possible, that\n*  abs(RT1R) >= abs(RT2R).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rt1r, rt1i, rt2r, rt2i, cs, sn, a, b, c, d = NumRu::Lapack.slanv2( a, b, c, d, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_a = argv[0];
  rblapack_b = argv[1];
  rblapack_c = argv[2];
  rblapack_d = argv[3];
  if (rb_options != Qnil) {
  }

  a = (real)NUM2DBL(rblapack_a);
  b = (real)NUM2DBL(rblapack_b);
  c = (real)NUM2DBL(rblapack_c);
  d = (real)NUM2DBL(rblapack_d);

  slanv2_(&a, &b, &c, &d, &rt1r, &rt1i, &rt2r, &rt2i, &cs, &sn);

  rblapack_rt1r = rb_float_new((double)rt1r);
  rblapack_rt1i = rb_float_new((double)rt1i);
  rblapack_rt2r = rb_float_new((double)rt2r);
  rblapack_rt2i = rb_float_new((double)rt2i);
  rblapack_cs = rb_float_new((double)cs);
  rblapack_sn = rb_float_new((double)sn);
  rblapack_a = rb_float_new((double)a);
  rblapack_b = rb_float_new((double)b);
  rblapack_c = rb_float_new((double)c);
  rblapack_d = rb_float_new((double)d);
  return rb_ary_new3(10, rblapack_rt1r, rblapack_rt1i, rblapack_rt2r, rblapack_rt2i, rblapack_cs, rblapack_sn, rblapack_a, rblapack_b, rblapack_c, rblapack_d);
}

void
init_lapack_slanv2(VALUE mLapack){
  rb_define_module_function(mLapack, "slanv2", rblapack_slanv2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
