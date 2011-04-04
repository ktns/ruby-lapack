#include "rb_lapack.h"

extern VOID csyequb_(char *uplo, integer *n, complex *a, integer *lda, real *s, real *scond, real *amax, complex *work, integer *info);

static VALUE
rb_csyequb(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
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
  complex *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, scond, amax, info = NumRu::Lapack.csyequb( uplo, a)\n    or\n  NumRu::Lapack.csyequb  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CSYEQUB computes row and column scalings intended to equilibrate a\n*  symmetric matrix A and reduce its condition number\n*  (with respect to the two-norm).  S contains the scale factors,\n*  S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with\n*  elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This\n*  choice of S puts the condition number of B within a factor N of the\n*  smallest possible condition number over all possible diagonal\n*  scalings.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**T;\n*          = 'L':  Lower triangular, form is A = L*D*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input) COMPLEX array, dimension (LDA,N)\n*          The N-by-N symmetric matrix whose scaling\n*          factors are to be computed.  Only the diagonal elements of A\n*          are referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  S       (output) REAL array, dimension (N)\n*          If INFO = 0, S contains the scale factors for A.\n*\n*  SCOND   (output) REAL\n*          If INFO = 0, S contains the ratio of the smallest S(i) to\n*          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too\n*          large nor too small, it is not worth scaling by S.\n*\n*  AMAX    (output) REAL\n*          Absolute value of largest matrix element.  If AMAX is very\n*          close to overflow or very close to underflow, the matrix\n*          should be scaled.\n*\n*  WORK    (workspace) COMPLEX array, dimension (3*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the i-th diagonal element is nonpositive.\n*\n\n*  Further Details\n*  ======= =======\n*\n*  Reference: Livne, O.E. and Golub, G.H., \"Scaling by Binormalization\",\n*  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004.\n*  DOI 10.1023/B:NUMA.0000016606.32820.69\n*  Tech report version: http://ruready.utah.edu/archive/papers/bin.pdf\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  work = ALLOC_N(complex, (3*n));

  csyequb_(&uplo, &n, a, &lda, s, &scond, &amax, work, &info);

  free(work);
  rb_scond = rb_float_new((double)scond);
  rb_amax = rb_float_new((double)amax);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_s, rb_scond, rb_amax, rb_info);
}

void
init_lapack_csyequb(VALUE mLapack){
  rb_define_module_function(mLapack, "csyequb", rb_csyequb, -1);
}
