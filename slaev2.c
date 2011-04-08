#include "rb_lapack.h"

extern VOID slaev2_(real *a, real *b, real *c, real *rt1, real *rt2, real *cs1, real *sn1);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaev2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  real a; 
  VALUE rblapack_b;
  real b; 
  VALUE rblapack_c;
  real c; 
  VALUE rblapack_rt1;
  real rt1; 
  VALUE rblapack_rt2;
  real rt2; 
  VALUE rblapack_cs1;
  real cs1; 
  VALUE rblapack_sn1;
  real sn1; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rt1, rt2, cs1, sn1 = NumRu::Lapack.slaev2( a, b, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAEV2( A, B, C, RT1, RT2, CS1, SN1 )\n\n*  Purpose\n*  =======\n*\n*  SLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix\n*     [  A   B  ]\n*     [  B   C  ].\n*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the\n*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right\n*  eigenvector for RT1, giving the decomposition\n*\n*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]\n*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input) REAL\n*          The (1,1) element of the 2-by-2 matrix.\n*\n*  B       (input) REAL\n*          The (1,2) element and the conjugate of the (2,1) element of\n*          the 2-by-2 matrix.\n*\n*  C       (input) REAL\n*          The (2,2) element of the 2-by-2 matrix.\n*\n*  RT1     (output) REAL\n*          The eigenvalue of larger absolute value.\n*\n*  RT2     (output) REAL\n*          The eigenvalue of smaller absolute value.\n*\n*  CS1     (output) REAL\n*  SN1     (output) REAL\n*          The vector (CS1, SN1) is a unit right eigenvector for RT1.\n*\n\n*  Further Details\n*  ===============\n*\n*  RT1 is accurate to a few ulps barring over/underflow.\n*\n*  RT2 may be inaccurate if there is massive cancellation in the\n*  determinant A*C-B*B; higher precision or correctly rounded or\n*  correctly truncated arithmetic would be needed to compute RT2\n*  accurately in all cases.\n*\n*  CS1 and SN1 are accurate to a few ulps barring over/underflow.\n*\n*  Overflow is possible only if RT1 is within a factor of 5 of overflow.\n*  Underflow is harmless if the input data is 0 or exceeds\n*     underflow_threshold / macheps.\n*\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rt1, rt2, cs1, sn1 = NumRu::Lapack.slaev2( a, b, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_a = argv[0];
  rblapack_b = argv[1];
  rblapack_c = argv[2];
  if (rb_options != Qnil) {
  }

  a = (real)NUM2DBL(rblapack_a);
  b = (real)NUM2DBL(rblapack_b);
  c = (real)NUM2DBL(rblapack_c);

  slaev2_(&a, &b, &c, &rt1, &rt2, &cs1, &sn1);

  rblapack_rt1 = rb_float_new((double)rt1);
  rblapack_rt2 = rb_float_new((double)rt2);
  rblapack_cs1 = rb_float_new((double)cs1);
  rblapack_sn1 = rb_float_new((double)sn1);
  return rb_ary_new3(4, rblapack_rt1, rblapack_rt2, rblapack_cs1, rblapack_sn1);
}

void
init_lapack_slaev2(VALUE mLapack){
  rb_define_module_function(mLapack, "slaev2", rblapack_slaev2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
