#include "rb_lapack.h"

static VALUE
rb_cpbequ(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_s;
  real *s; 
  VALUE rb_scond;
  real scond; 
  VALUE rb_amax;
  real amax; 
  VALUE rb_info;
  integer info; 

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, scond, amax, info = NumRu::Lapack.cpbequ( uplo, kd, ab)\n    or\n  NumRu::Lapack.cpbequ  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )\n\n*  Purpose\n*  =======\n*\n*  CPBEQU computes row and column scalings intended to equilibrate a\n*  Hermitian positive definite band matrix A and reduce its condition\n*  number (with respect to the two-norm).  S contains the scale factors,\n*  S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with\n*  elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This\n*  choice of S puts the condition number of B within a factor N of the\n*  smallest possible condition number over all possible diagonal\n*  scalings.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangular of A is stored;\n*          = 'L':  Lower triangular of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input) COMPLEX array, dimension (LDAB,N)\n*          The upper or lower triangle of the Hermitian band matrix A,\n*          stored in the first KD+1 rows of the array.  The j-th column\n*          of A is stored in the j-th column of the array AB as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*\n*  LDAB     (input) INTEGER\n*          The leading dimension of the array A.  LDAB >= KD+1.\n*\n*  S       (output) REAL array, dimension (N)\n*          If INFO = 0, S contains the scale factors for A.\n*\n*  SCOND   (output) REAL\n*          If INFO = 0, S contains the ratio of the smallest S(i) to\n*          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too\n*          large nor too small, it is not worth scaling by S.\n*\n*  AMAX    (output) REAL\n*          Absolute value of largest matrix element.  If AMAX is very\n*          close to overflow or very close to underflow, the matrix\n*          should be scaled.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, the i-th diagonal element is nonpositive.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_kd = argv[1];
  rb_ab = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  kd = NUM2INT(rb_kd);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);

  cpbequ_(&uplo, &n, &kd, ab, &ldab, s, &scond, &amax, &info);

  rb_scond = rb_float_new((double)scond);
  rb_amax = rb_float_new((double)amax);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_s, rb_scond, rb_amax, rb_info);
}

void
init_lapack_cpbequ(VALUE mLapack){
  rb_define_module_function(mLapack, "cpbequ", rb_cpbequ, -1);
}
