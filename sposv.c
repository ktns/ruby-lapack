#include "rb_lapack.h"

extern VOID sposv_(char *uplo, integer *n, integer *nrhs, real *a, integer *lda, real *b, integer *ldb, integer *info);

static VALUE
rb_sposv(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a, b = NumRu::Lapack.sposv( uplo, a, b)\n    or\n  NumRu::Lapack.sposv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPOSV computes the solution to a real system of linear equations\n*     A * X = B,\n*  where A is an N-by-N symmetric positive definite matrix and X and B\n*  are N-by-NRHS matrices.\n*\n*  The Cholesky decomposition is used to factor A as\n*     A = U**T* U,  if UPLO = 'U', or\n*     A = L * L**T,  if UPLO = 'L',\n*  where U is an upper triangular matrix and L is a lower triangular\n*  matrix.  The factored form of A is then used to solve the system of\n*  equations A * X = B.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading\n*          N-by-N upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading N-by-N lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*          On exit, if INFO = 0, the factor U or L from the Cholesky\n*          factorization A = U**T*U or A = L*L**T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the leading minor of order i of A is not\n*                positive definite, so the factorization could not be\n*                completed, and the solution has not been computed.\n*\n\n*  =====================================================================\n*\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SPOTRF, SPOTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
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

  sposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_a, rb_b);
}

void
init_lapack_sposv(VALUE mLapack){
  rb_define_module_function(mLapack, "sposv", rb_sposv, -1);
}
