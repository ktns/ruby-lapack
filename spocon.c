#include "rb_lapack.h"

static VALUE
rb_spocon(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_anorm;
  real anorm; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_info;
  integer info; 
  real *work;
  integer *iwork;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.spocon( uplo, a, anorm)\n    or\n  NumRu::Lapack.spocon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPOCON estimates the reciprocal of the condition number (in the \n*  1-norm) of a real symmetric positive definite matrix using the\n*  Cholesky factorization A = U**T*U or A = L*L**T computed by SPOTRF.\n*\n*  An estimate is obtained for norm(inv(A)), and the reciprocal of the\n*  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input) REAL array, dimension (LDA,N)\n*          The triangular factor U or L from the Cholesky factorization\n*          A = U**T*U or A = L*L**T, as computed by SPOTRF.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  ANORM   (input) REAL\n*          The 1-norm (or infinity-norm) of the symmetric matrix A.\n*\n*  RCOND   (output) REAL\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an\n*          estimate of the 1-norm of inv(A) computed in this routine.\n*\n*  WORK    (workspace) REAL array, dimension (3*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_anorm = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  anorm = (real)NUM2DBL(rb_anorm);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  work = ALLOC_N(real, (3*n));
  iwork = ALLOC_N(integer, (n));

  spocon_(&uplo, &n, a, &lda, &anorm, &rcond, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_spocon(VALUE mLapack){
  rb_define_module_function(mLapack, "spocon", rb_spocon, -1);
}
