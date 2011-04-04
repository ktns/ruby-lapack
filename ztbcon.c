#include "rb_lapack.h"

extern VOID ztbcon_(char *norm, char *uplo, char *diag, integer *n, integer *kd, doublecomplex *ab, integer *ldab, doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE
rb_ztbcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_info;
  integer info; 
  doublecomplex *work;
  doublereal *rwork;

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.ztbcon( norm, uplo, diag, kd, ab)\n    or\n  NumRu::Lapack.ztbcon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTBCON estimates the reciprocal of the condition number of a\n*  triangular band matrix A, in either the 1-norm or the infinity-norm.\n*\n*  The norm of A is computed and an estimate is obtained for\n*  norm(inv(A)), then the reciprocal of the condition number is\n*  computed as\n*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies whether the 1-norm condition number or the\n*          infinity-norm condition number is required:\n*          = '1' or 'O':  1-norm;\n*          = 'I':         Infinity-norm.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals or subdiagonals of the\n*          triangular band matrix A.  KD >= 0.\n*\n*  AB      (input) COMPLEX*16 array, dimension (LDAB,N)\n*          The upper or lower triangular band matrix A, stored in the\n*          first kd+1 rows of the array. The j-th column of A is stored\n*          in the j-th column of the array AB as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*          If DIAG = 'U', the diagonal elements of A are not referenced\n*          and are assumed to be 1.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD+1.\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(norm(A) * norm(inv(A))).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (2*N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_kd = argv[3];
  rb_ab = argv[4];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  diag = StringValueCStr(rb_diag)[0];
  kd = NUM2INT(rb_kd);
  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (n));

  ztbcon_(&norm, &uplo, &diag, &n, &kd, ab, &ldab, &rcond, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_ztbcon(VALUE mLapack){
  rb_define_module_function(mLapack, "ztbcon", rb_ztbcon, -1);
}
