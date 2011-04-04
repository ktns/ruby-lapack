#include "rb_lapack.h"

extern VOID slanv2_(real *a, real *b, real *c, real *d, real *rt1r, real *rt1i, real *rt2r, real *rt2i, real *cs, real *sn);

static VALUE
rb_slanv2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real a; 
  VALUE rb_b;
  real b; 
  VALUE rb_c;
  real c; 
  VALUE rb_d;
  real d; 
  VALUE rb_rt1r;
  real rt1r; 
  VALUE rb_rt1i;
  real rt1i; 
  VALUE rb_rt2r;
  real rt2r; 
  VALUE rb_rt2i;
  real rt2i; 
  VALUE rb_cs;
  real cs; 
  VALUE rb_sn;
  real sn; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rt1r, rt1i, rt2r, rt2i, cs, sn, a, b, c, d = NumRu::Lapack.slanv2( a, b, c, d)\n    or\n  NumRu::Lapack.slanv2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )\n\n*  Purpose\n*  =======\n*\n*  SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric\n*  matrix in standard form:\n*\n*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]\n*       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]\n*\n*  where either\n*  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or\n*  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex\n*  conjugate eigenvalues.\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input/output) REAL            \n*  B       (input/output) REAL            \n*  C       (input/output) REAL            \n*  D       (input/output) REAL            \n*          On entry, the elements of the input matrix.\n*          On exit, they are overwritten by the elements of the\n*          standardised Schur form.\n*\n*  RT1R    (output) REAL \n*  RT1I    (output) REAL            \n*  RT2R    (output) REAL            \n*  RT2I    (output) REAL            \n*          The real and imaginary parts of the eigenvalues. If the\n*          eigenvalues are a complex conjugate pair, RT1I > 0.\n*\n*  CS      (output) REAL            \n*  SN      (output) REAL            \n*          Parameters of the rotation matrix.\n*\n\n*  Further Details\n*  ===============\n*\n*  Modified by V. Sima, Research Institute for Informatics, Bucharest,\n*  Romania, to reduce the risk of cancellation errors,\n*  when computing real eigenvalues, and to ensure, if possible, that\n*  abs(RT1R) >= abs(RT2R).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_c = argv[2];
  rb_d = argv[3];

  a = (real)NUM2DBL(rb_a);
  b = (real)NUM2DBL(rb_b);
  c = (real)NUM2DBL(rb_c);
  d = (real)NUM2DBL(rb_d);

  slanv2_(&a, &b, &c, &d, &rt1r, &rt1i, &rt2r, &rt2i, &cs, &sn);

  rb_rt1r = rb_float_new((double)rt1r);
  rb_rt1i = rb_float_new((double)rt1i);
  rb_rt2r = rb_float_new((double)rt2r);
  rb_rt2i = rb_float_new((double)rt2i);
  rb_cs = rb_float_new((double)cs);
  rb_sn = rb_float_new((double)sn);
  rb_a = rb_float_new((double)a);
  rb_b = rb_float_new((double)b);
  rb_c = rb_float_new((double)c);
  rb_d = rb_float_new((double)d);
  return rb_ary_new3(10, rb_rt1r, rb_rt1i, rb_rt2r, rb_rt2i, rb_cs, rb_sn, rb_a, rb_b, rb_c, rb_d);
}

void
init_lapack_slanv2(VALUE mLapack){
  rb_define_module_function(mLapack, "slanv2", rb_slanv2, -1);
}
