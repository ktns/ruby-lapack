#include "rb_lapack.h"

extern VOID dsycon_(char *uplo, integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dsycon(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_anorm;
  doublereal anorm; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_info;
  integer info; 
  doublereal *work;
  integer *iwork;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.dsycon( uplo, a, ipiv, anorm)\n    or\n  NumRu::Lapack.dsycon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSYCON estimates the reciprocal of the condition number (in the\n*  1-norm) of a real symmetric matrix A using the factorization\n*  A = U*D*U**T or A = L*D*L**T computed by DSYTRF.\n*\n*  An estimate is obtained for norm(inv(A)), and the reciprocal of the\n*  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**T;\n*          = 'L':  Lower triangular, form is A = L*D*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*          The block diagonal matrix D and the multipliers used to\n*          obtain the factor U or L as computed by DSYTRF.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D\n*          as determined by DSYTRF.\n*\n*  ANORM   (input) DOUBLE PRECISION\n*          The 1-norm of the original matrix A.\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an\n*          estimate of the 1-norm of inv(A) computed in this routine.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_ipiv = argv[2];
  rb_anorm = argv[3];

  anorm = NUM2DBL(rb_anorm);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (3th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  uplo = StringValueCStr(rb_uplo)[0];
  work = ALLOC_N(doublereal, (2*n));
  iwork = ALLOC_N(integer, (n));

  dsycon_(&uplo, &n, a, &lda, ipiv, &anorm, &rcond, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_dsycon(VALUE mLapack){
  rb_define_module_function(mLapack, "dsycon", rb_dsycon, -1);
}
