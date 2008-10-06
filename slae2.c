#include "rb_lapack.h"

static VALUE
rb_slae2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real a; 
  VALUE rb_b;
  real b; 
  VALUE rb_c;
  real c; 
  VALUE rb_rt1;
  real rt1; 
  VALUE rb_rt2;
  real rt2; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rt1, rt2 = NumRu::Lapack.slae2( a, b, c)\n    or\n  NumRu::Lapack.slae2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAE2( A, B, C, RT1, RT2 )\n\n*  Purpose\n*  =======\n*\n*  SLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix\n*     [  A   B  ]\n*     [  B   C  ].\n*  On return, RT1 is the eigenvalue of larger absolute value, and RT2\n*  is the eigenvalue of smaller absolute value.\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input) REAL\n*          The (1,1) element of the 2-by-2 matrix.\n*\n*  B       (input) REAL\n*          The (1,2) and (2,1) elements of the 2-by-2 matrix.\n*\n*  C       (input) REAL\n*          The (2,2) element of the 2-by-2 matrix.\n*\n*  RT1     (output) REAL\n*          The eigenvalue of larger absolute value.\n*\n*  RT2     (output) REAL\n*          The eigenvalue of smaller absolute value.\n*\n\n*  Further Details\n*  ===============\n*\n*  RT1 is accurate to a few ulps barring over/underflow.\n*\n*  RT2 may be inaccurate if there is massive cancellation in the\n*  determinant A*C-B*B; higher precision or correctly rounded or\n*  correctly truncated arithmetic would be needed to compute RT2\n*  accurately in all cases.\n*\n*  Overflow is possible only if RT1 is within a factor of 5 of overflow.\n*  Underflow is harmless if the input data is 0 or exceeds\n*     underflow_threshold / macheps.\n*\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_c = argv[2];

  a = (real)NUM2DBL(rb_a);
  b = (real)NUM2DBL(rb_b);
  c = (real)NUM2DBL(rb_c);

  slae2_(&a, &b, &c, &rt1, &rt2);

  rb_rt1 = rb_float_new((double)rt1);
  rb_rt2 = rb_float_new((double)rt2);
  return rb_ary_new3(2, rb_rt1, rb_rt2);
}

void
init_lapack_slae2(VALUE mLapack){
  rb_define_module_function(mLapack, "slae2", rb_slae2, -1);
}
