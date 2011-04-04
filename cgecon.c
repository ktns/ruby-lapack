#include "rb_lapack.h"

extern VOID cgecon_(char *norm, integer *n, complex *a, integer *lda, real *anorm, real *rcond, complex *work, real *rwork, integer *info);

static VALUE
rb_cgecon(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_anorm;
  real anorm; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_info;
  integer info; 
  complex *work;
  real *rwork;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.cgecon( norm, a, anorm)\n    or\n  NumRu::Lapack.cgecon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGECON estimates the reciprocal of the condition number of a general\n*  complex matrix A, in either the 1-norm or the infinity-norm, using\n*  the LU factorization computed by CGETRF.\n*\n*  An estimate is obtained for norm(inv(A)), and the reciprocal of the\n*  condition number is computed as\n*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies whether the 1-norm condition number or the\n*          infinity-norm condition number is required:\n*          = '1' or 'O':  1-norm;\n*          = 'I':         Infinity-norm.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input) COMPLEX array, dimension (LDA,N)\n*          The factors L and U from the factorization A = P*L*U\n*          as computed by CGETRF.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  ANORM   (input) REAL\n*          If NORM = '1' or 'O', the 1-norm of the original matrix A.\n*          If NORM = 'I', the infinity-norm of the original matrix A.\n*\n*  RCOND   (output) REAL\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(norm(A) * norm(inv(A))).\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  RWORK   (workspace) REAL array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_norm = argv[0];
  rb_a = argv[1];
  rb_anorm = argv[2];

  anorm = (real)NUM2DBL(rb_anorm);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  norm = StringValueCStr(rb_norm)[0];
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (2*n));

  cgecon_(&norm, &n, a, &lda, &anorm, &rcond, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_cgecon(VALUE mLapack){
  rb_define_module_function(mLapack, "cgecon", rb_cgecon, -1);
}
