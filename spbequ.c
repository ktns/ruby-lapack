#include "rb_lapack.h"

extern VOID spbequ_(char *uplo, integer *n, integer *kd, real *ab, integer *ldab, real *s, real *scond, real *amax, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_spbequ(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_kd;
  integer kd; 
  VALUE rblapack_ab;
  real *ab; 
  VALUE rblapack_s;
  real *s; 
  VALUE rblapack_scond;
  real scond; 
  VALUE rblapack_amax;
  real amax; 
  VALUE rblapack_info;
  integer info; 

  integer ldab;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, scond, amax, info = NumRu::Lapack.spbequ( uplo, kd, ab, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPBEQU computes row and column scalings intended to equilibrate a\n*  symmetric positive definite band matrix A and reduce its condition\n*  number (with respect to the two-norm).  S contains the scale factors,\n*  S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with\n*  elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This\n*  choice of S puts the condition number of B within a factor N of the\n*  smallest possible condition number over all possible diagonal\n*  scalings.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangular of A is stored;\n*          = 'L':  Lower triangular of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input) REAL array, dimension (LDAB,N)\n*          The upper or lower triangle of the symmetric band matrix A,\n*          stored in the first KD+1 rows of the array.  The j-th column\n*          of A is stored in the j-th column of the array AB as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*\n*  LDAB     (input) INTEGER\n*          The leading dimension of the array A.  LDAB >= KD+1.\n*\n*  S       (output) REAL array, dimension (N)\n*          If INFO = 0, S contains the scale factors for A.\n*\n*  SCOND   (output) REAL\n*          If INFO = 0, S contains the ratio of the smallest S(i) to\n*          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too\n*          large nor too small, it is not worth scaling by S.\n*\n*  AMAX    (output) REAL\n*          Absolute value of largest matrix element.  If AMAX is very\n*          close to overflow or very close to underflow, the matrix\n*          should be scaled.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, the i-th diagonal element is nonpositive.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, scond, amax, info = NumRu::Lapack.spbequ( uplo, kd, ab, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_uplo = argv[0];
  rblapack_kd = argv[1];
  rblapack_ab = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_ab);
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SFLOAT)
    rblapack_ab = na_change_type(rblapack_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rblapack_ab, real*);
  kd = NUM2INT(rblapack_kd);
  uplo = StringValueCStr(rblapack_uplo)[0];
  {
    int shape[1];
    shape[0] = n;
    rblapack_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rblapack_s, real*);

  spbequ_(&uplo, &n, &kd, ab, &ldab, s, &scond, &amax, &info);

  rblapack_scond = rb_float_new((double)scond);
  rblapack_amax = rb_float_new((double)amax);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_s, rblapack_scond, rblapack_amax, rblapack_info);
}

void
init_lapack_spbequ(VALUE mLapack){
  rb_define_module_function(mLapack, "spbequ", rblapack_spbequ, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
