#include "rb_lapack.h"

extern VOID spbsv_(char *uplo, integer *n, integer *kd, integer *nrhs, real *ab, integer *ldab, real *b, integer *ldb, integer *info);

static VALUE
rb_spbsv(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_b;
  real *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  real *ab_out__;
  VALUE rb_b_out__;
  real *b_out__;

  integer ldab;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, ab, b = NumRu::Lapack.spbsv( uplo, kd, ab, b)\n    or\n  NumRu::Lapack.spbsv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPBSV computes the solution to a real system of linear equations\n*     A * X = B,\n*  where A is an N-by-N symmetric positive definite band matrix and X\n*  and B are N-by-NRHS matrices.\n*\n*  The Cholesky decomposition is used to factor A as\n*     A = U**T * U,  if UPLO = 'U', or\n*     A = L * L**T,  if UPLO = 'L',\n*  where U is an upper triangular band matrix, and L is a lower\n*  triangular band matrix, with the same number of superdiagonals or\n*  subdiagonals as A.  The factored form of A is then used to solve the\n*  system of equations A * X = B.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  AB      (input/output) REAL array, dimension (LDAB,N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first KD+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD).\n*          See below for further details.\n*\n*          On exit, if INFO = 0, the triangular factor U or L from the\n*          Cholesky factorization A = U**T*U or A = L*L**T of the band\n*          matrix A, in the same storage format as A.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD+1.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the leading minor of order i of A is not\n*                positive definite, so the factorization could not be\n*                completed, and the solution has not been computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  The band storage scheme is illustrated by the following example, when\n*  N = 6, KD = 2, and UPLO = 'U':\n*\n*  On entry:                       On exit:\n*\n*      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46\n*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56\n*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66\n*\n*  Similarly, if UPLO = 'L' the format of A is as follows:\n*\n*  On entry:                       On exit:\n*\n*     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66\n*     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *\n*     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *\n*\n*  Array elements marked * are not used by the routine.\n*\n*  =====================================================================\n*\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SPBTRF, SPBTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_kd = argv[1];
  rb_ab = argv[2];
  rb_b = argv[3];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  kd = NUM2INT(rb_kd);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  spbsv_(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_ab, rb_b);
}

void
init_lapack_spbsv(VALUE mLapack){
  rb_define_module_function(mLapack, "spbsv", rb_spbsv, -1);
}
