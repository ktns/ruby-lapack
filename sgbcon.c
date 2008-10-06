#include "rb_lapack.h"

static VALUE
rb_sgbcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_anorm;
  real anorm; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_info;
  integer info; 
  real *work;
  integer *iwork;

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.sgbcon( norm, kl, ku, ab, ipiv, anorm)\n    or\n  NumRu::Lapack.sgbcon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGBCON estimates the reciprocal of the condition number of a real\n*  general band matrix A, in either the 1-norm or the infinity-norm,\n*  using the LU factorization computed by SGBTRF.\n*\n*  An estimate is obtained for norm(inv(A)), and the reciprocal of the\n*  condition number is computed as\n*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies whether the 1-norm condition number or the\n*          infinity-norm condition number is required:\n*          = '1' or 'O':  1-norm;\n*          = 'I':         Infinity-norm.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  AB      (input) REAL array, dimension (LDAB,N)\n*          Details of the LU factorization of the band matrix A, as\n*          computed by SGBTRF.  U is stored as an upper triangular band\n*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and\n*          the multipliers used during the factorization are stored in\n*          rows KL+KU+2 to 2*KL+KU+1.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= N, row i of the matrix was\n*          interchanged with row IPIV(i).\n*\n*  ANORM   (input) REAL\n*          If NORM = '1' or 'O', the 1-norm of the original matrix A.\n*          If NORM = 'I', the infinity-norm of the original matrix A.\n*\n*  RCOND   (output) REAL\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(norm(A) * norm(inv(A))).\n*\n*  WORK    (workspace) REAL array, dimension (3*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_norm = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];
  rb_ipiv = argv[4];
  rb_anorm = argv[5];

  norm = StringValueCStr(rb_norm)[0];
  kl = NUM2INT(rb_kl);
  ku = NUM2INT(rb_ku);
  anorm = (real)NUM2DBL(rb_anorm);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ipiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be the same as shape 1 of ab");
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  work = ALLOC_N(real, (3*n));
  iwork = ALLOC_N(integer, (n));

  sgbcon_(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, &rcond, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_sgbcon(VALUE mLapack){
  rb_define_module_function(mLapack, "sgbcon", rb_sgbcon, -1);
}
