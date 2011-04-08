#include "rb_lapack.h"

extern VOID dlaev2_(doublereal *a, doublereal *b, doublereal *c, doublereal *rt1, doublereal *rt2, doublereal *cs1, doublereal *sn1);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlaev2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  doublereal a; 
  VALUE rblapack_b;
  doublereal b; 
  VALUE rblapack_c;
  doublereal c; 
  VALUE rblapack_rt1;
  doublereal rt1; 
  VALUE rblapack_rt2;
  doublereal rt2; 
  VALUE rblapack_cs1;
  doublereal cs1; 
  VALUE rblapack_sn1;
  doublereal sn1; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rt1, rt2, cs1, sn1 = NumRu::Lapack.dlaev2( a, b, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )\n\n*  Purpose\n*  =======\n*\n*  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix\n*     [  A   B  ]\n*     [  B   C  ].\n*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the\n*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right\n*  eigenvector for RT1, giving the decomposition\n*\n*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]\n*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input) DOUBLE PRECISION\n*          The (1,1) element of the 2-by-2 matrix.\n*\n*  B       (input) DOUBLE PRECISION\n*          The (1,2) element and the conjugate of the (2,1) element of\n*          the 2-by-2 matrix.\n*\n*  C       (input) DOUBLE PRECISION\n*          The (2,2) element of the 2-by-2 matrix.\n*\n*  RT1     (output) DOUBLE PRECISION\n*          The eigenvalue of larger absolute value.\n*\n*  RT2     (output) DOUBLE PRECISION\n*          The eigenvalue of smaller absolute value.\n*\n*  CS1     (output) DOUBLE PRECISION\n*  SN1     (output) DOUBLE PRECISION\n*          The vector (CS1, SN1) is a unit right eigenvector for RT1.\n*\n\n*  Further Details\n*  ===============\n*\n*  RT1 is accurate to a few ulps barring over/underflow.\n*\n*  RT2 may be inaccurate if there is massive cancellation in the\n*  determinant A*C-B*B; higher precision or correctly rounded or\n*  correctly truncated arithmetic would be needed to compute RT2\n*  accurately in all cases.\n*\n*  CS1 and SN1 are accurate to a few ulps barring over/underflow.\n*\n*  Overflow is possible only if RT1 is within a factor of 5 of overflow.\n*  Underflow is harmless if the input data is 0 or exceeds\n*     underflow_threshold / macheps.\n*\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rt1, rt2, cs1, sn1 = NumRu::Lapack.dlaev2( a, b, c, [:usage => usage, :help => help])\n");
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

  a = NUM2DBL(rblapack_a);
  b = NUM2DBL(rblapack_b);
  c = NUM2DBL(rblapack_c);

  dlaev2_(&a, &b, &c, &rt1, &rt2, &cs1, &sn1);

  rblapack_rt1 = rb_float_new((double)rt1);
  rblapack_rt2 = rb_float_new((double)rt2);
  rblapack_cs1 = rb_float_new((double)cs1);
  rblapack_sn1 = rb_float_new((double)sn1);
  return rb_ary_new3(4, rblapack_rt1, rblapack_rt2, rblapack_cs1, rblapack_sn1);
}

void
init_lapack_dlaev2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaev2", rblapack_dlaev2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
