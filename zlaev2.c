#include "rb_lapack.h"

static VALUE
rb_zlaev2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublecomplex a; 
  VALUE rb_b;
  doublecomplex b; 
  VALUE rb_c;
  doublecomplex c; 
  VALUE rb_rt1;
  doublereal rt1; 
  VALUE rb_rt2;
  doublereal rt2; 
  VALUE rb_cs1;
  doublereal cs1; 
  VALUE rb_sn1;
  doublecomplex sn1; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rt1, rt2, cs1, sn1 = NumRu::Lapack.zlaev2( a, b, c)\n    or\n  NumRu::Lapack.zlaev2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAEV2( A, B, C, RT1, RT2, CS1, SN1 )\n\n*  Purpose\n*  =======\n*\n*  ZLAEV2 computes the eigendecomposition of a 2-by-2 Hermitian matrix\n*     [  A         B  ]\n*     [  CONJG(B)  C  ].\n*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the\n*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right\n*  eigenvector for RT1, giving the decomposition\n*\n*  [ CS1  CONJG(SN1) ] [    A     B ] [ CS1 -CONJG(SN1) ] = [ RT1  0  ]\n*  [-SN1     CS1     ] [ CONJG(B) C ] [ SN1     CS1     ]   [  0  RT2 ].\n*\n\n*  Arguments\n*  =========\n*\n*  A      (input) COMPLEX*16\n*         The (1,1) element of the 2-by-2 matrix.\n*\n*  B      (input) COMPLEX*16\n*         The (1,2) element and the conjugate of the (2,1) element of\n*         the 2-by-2 matrix.\n*\n*  C      (input) COMPLEX*16\n*         The (2,2) element of the 2-by-2 matrix.\n*\n*  RT1    (output) DOUBLE PRECISION\n*         The eigenvalue of larger absolute value.\n*\n*  RT2    (output) DOUBLE PRECISION\n*         The eigenvalue of smaller absolute value.\n*\n*  CS1    (output) DOUBLE PRECISION\n*  SN1    (output) COMPLEX*16\n*         The vector (CS1, SN1) is a unit right eigenvector for RT1.\n*\n\n*  Further Details\n*  ===============\n*\n*  RT1 is accurate to a few ulps barring over/underflow.\n*\n*  RT2 may be inaccurate if there is massive cancellation in the\n*  determinant A*C-B*B; higher precision or correctly rounded or\n*  correctly truncated arithmetic would be needed to compute RT2\n*  accurately in all cases.\n*\n*  CS1 and SN1 are accurate to a few ulps barring over/underflow.\n*\n*  Overflow is possible only if RT1 is within a factor of 5 of overflow.\n*  Underflow is harmless if the input data is 0 or exceeds\n*     underflow_threshold / macheps.\n*\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_c = argv[2];

  a.r = NUM2DBL(rb_funcall(rb_a, rb_intern("real"), 0));
  a.i = NUM2DBL(rb_funcall(rb_a, rb_intern("imag"), 0));
  b.r = NUM2DBL(rb_funcall(rb_b, rb_intern("real"), 0));
  b.i = NUM2DBL(rb_funcall(rb_b, rb_intern("imag"), 0));
  c.r = NUM2DBL(rb_funcall(rb_c, rb_intern("real"), 0));
  c.i = NUM2DBL(rb_funcall(rb_c, rb_intern("imag"), 0));

  zlaev2_(&a, &b, &c, &rt1, &rt2, &cs1, &sn1);

  rb_rt1 = rb_float_new((double)rt1);
  rb_rt2 = rb_float_new((double)rt2);
  rb_cs1 = rb_float_new((double)cs1);
  rb_sn1 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(sn1.r)), rb_float_new((double)(sn1.i)));
  return rb_ary_new3(4, rb_rt1, rb_rt2, rb_cs1, rb_sn1);
}

void
init_lapack_zlaev2(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaev2", rb_zlaev2, -1);
}
