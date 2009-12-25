#include "rb_lapack.h"

static VALUE
rb_zsyequb(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_work;
  doublereal work; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_scond;
  doublereal scond; 
  VALUE rb_amax;
  doublereal amax; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, scond, amax, info = NumRu::Lapack.zsyequb( uplo, a, work)\n    or\n  NumRu::Lapack.zsyequb  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZSYEQUB computes row and column scalings intended to equilibrate a\n*  symmetric matrix A and reduce its condition number\n*  (with respect to the two-norm).  S contains the scale factors,\n*  S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with\n*  elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This\n*  choice of S puts the condition number of B within a factor N of the\n*  smallest possible condition number over all possible diagonal\n*  scalings.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA,N)\n*          The N-by-N symmetric matrix whose scaling\n*          factors are to be computed.  Only the diagonal elements of A\n*          are referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  S       (output) DOUBLE PRECISION array, dimension (N)\n*          If INFO = 0, S contains the scale factors for A.\n*\n*  SCOND   (output) DOUBLE PRECISION\n*          If INFO = 0, S contains the ratio of the smallest S(i) to\n*          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too\n*          large nor too small, it is not worth scaling by S.\n*\n*  AMAX    (output) DOUBLE PRECISION\n*          Absolute value of largest matrix element.  If AMAX is very\n*          close to overflow or very close to underflow, the matrix\n*          should be scaled.\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the i-th diagonal element is nonpositive.\n*\n\n*  Further Details\n*  ======= =======\n*\n*  Reference: Livne, O.E. and Golub, G.H., \"Scaling by Binormalization\",\n*  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004.\n*  DOI 10.1023/B:NUMA.0000016606.32820.69\n*  Tech report version: http://ruready.utah.edu/archive/papers/bin.pdf\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_work = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  work = NUM2DBL(rb_work);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, doublereal*);

  zsyequb_(&uplo, &n, a, &lda, s, &scond, &amax, &work, &info);

  rb_scond = rb_float_new((double)scond);
  rb_amax = rb_float_new((double)amax);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_s, rb_scond, rb_amax, rb_info);
}

void
init_lapack_zsyequb(VALUE mLapack){
  rb_define_module_function(mLapack, "zsyequb", rb_zsyequb, -1);
}
