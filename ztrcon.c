#include "rb_lapack.h"

extern VOID ztrcon_(char *norm, char *uplo, char *diag, integer *n, doublecomplex *a, integer *lda, doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE
rb_ztrcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_info;
  integer info; 
  doublecomplex *work;
  doublereal *rwork;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.ztrcon( norm, uplo, diag, a)\n    or\n  NumRu::Lapack.ztrcon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTRCON estimates the reciprocal of the condition number of a\n*  triangular matrix A, in either the 1-norm or the infinity-norm.\n*\n*  The norm of A is computed and an estimate is obtained for\n*  norm(inv(A)), then the reciprocal of the condition number is\n*  computed as\n*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies whether the 1-norm condition number or the\n*          infinity-norm condition number is required:\n*          = '1' or 'O':  1-norm;\n*          = 'I':         Infinity-norm.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA,N)\n*          The triangular matrix A.  If UPLO = 'U', the leading N-by-N\n*          upper triangular part of the array A contains the upper\n*          triangular matrix, and the strictly lower triangular part of\n*          A is not referenced.  If UPLO = 'L', the leading N-by-N lower\n*          triangular part of the array A contains the lower triangular\n*          matrix, and the strictly upper triangular part of A is not\n*          referenced.  If DIAG = 'U', the diagonal elements of A are\n*          also not referenced and are assumed to be 1.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(norm(A) * norm(inv(A))).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (2*N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_a = argv[3];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  diag = StringValueCStr(rb_diag)[0];
  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (n));

  ztrcon_(&norm, &uplo, &diag, &n, a, &lda, &rcond, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_ztrcon(VALUE mLapack){
  rb_define_module_function(mLapack, "ztrcon", rb_ztrcon, -1);
}
