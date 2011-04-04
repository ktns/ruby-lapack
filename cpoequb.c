#include "rb_lapack.h"

extern VOID cpoequb_(integer *n, complex *a, integer *lda, real *s, real *scond, real *amax, integer *info);

static VALUE
rb_cpoequb(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  complex *a; 
  VALUE rb_s;
  real *s; 
  VALUE rb_scond;
  real scond; 
  VALUE rb_amax;
  real amax; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, scond, amax, info = NumRu::Lapack.cpoequb( a)\n    or\n  NumRu::Lapack.cpoequb  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO )\n\n*  Purpose\n*  =======\n*\n*  CPOEQUB computes row and column scalings intended to equilibrate a\n*  symmetric positive definite matrix A and reduce its condition number\n*  (with respect to the two-norm).  S contains the scale factors,\n*  S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with\n*  elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This\n*  choice of S puts the condition number of B within a factor N of the\n*  smallest possible condition number over all possible diagonal\n*  scalings.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input) COMPLEX array, dimension (LDA,N)\n*          The N-by-N symmetric positive definite matrix whose scaling\n*          factors are to be computed.  Only the diagonal elements of A\n*          are referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  S       (output) REAL array, dimension (N)\n*          If INFO = 0, S contains the scale factors for A.\n*\n*  SCOND   (output) REAL\n*          If INFO = 0, S contains the ratio of the smallest S(i) to\n*          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too\n*          large nor too small, it is not worth scaling by S.\n*\n*  AMAX    (output) REAL\n*          Absolute value of largest matrix element.  If AMAX is very\n*          close to overflow or very close to underflow, the matrix\n*          should be scaled.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the i-th diagonal element is nonpositive.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_a = argv[0];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);

  cpoequb_(&n, a, &lda, s, &scond, &amax, &info);

  rb_scond = rb_float_new((double)scond);
  rb_amax = rb_float_new((double)amax);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_s, rb_scond, rb_amax, rb_info);
}

void
init_lapack_cpoequb(VALUE mLapack){
  rb_define_module_function(mLapack, "cpoequb", rb_cpoequb, -1);
}
