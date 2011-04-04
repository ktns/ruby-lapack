#include "rb_lapack.h"

extern VOID stpcon_(char *norm, char *uplo, char *diag, integer *n, real *ap, real *rcond, real *work, integer *iwork, integer *info);

static VALUE
rb_stpcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_info;
  integer info; 
  real *work;
  integer *iwork;

  integer ldap;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.stpcon( norm, uplo, diag, ap)\n    or\n  NumRu::Lapack.stpcon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE STPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  STPCON estimates the reciprocal of the condition number of a packed\n*  triangular matrix A, in either the 1-norm or the infinity-norm.\n*\n*  The norm of A is computed and an estimate is obtained for\n*  norm(inv(A)), then the reciprocal of the condition number is\n*  computed as\n*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies whether the 1-norm condition number or the\n*          infinity-norm condition number is required:\n*          = '1' or 'O':  1-norm;\n*          = 'I':         Infinity-norm.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input) REAL array, dimension (N*(N+1)/2)\n*          The upper or lower triangular matrix A, packed columnwise in\n*          a linear array.  The j-th column of A is stored in the array\n*          AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*          If DIAG = 'U', the diagonal elements of A are not referenced\n*          and are assumed to be 1.\n*\n*  RCOND   (output) REAL\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(norm(A) * norm(inv(A))).\n*\n*  WORK    (workspace) REAL array, dimension (3*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_ap = argv[3];

  diag = StringValueCStr(rb_diag)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  norm = StringValueCStr(rb_norm)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_SFLOAT)
    rb_ap = na_change_type(rb_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rb_ap, real*);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  work = ALLOC_N(real, (3*n));
  iwork = ALLOC_N(integer, (n));

  stpcon_(&norm, &uplo, &diag, &n, ap, &rcond, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_stpcon(VALUE mLapack){
  rb_define_module_function(mLapack, "stpcon", rb_stpcon, -1);
}
