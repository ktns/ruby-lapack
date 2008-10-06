#include "rb_lapack.h"

static VALUE
rb_sppcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_anorm;
  real anorm; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_info;
  integer info; 
  real *work;
  integer *iwork;

  integer ldap;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.sppcon( uplo, ap, anorm)\n    or\n  NumRu::Lapack.sppcon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPPCON( UPLO, N, AP, ANORM, RCOND, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPPCON estimates the reciprocal of the condition number (in the\n*  1-norm) of a real symmetric positive definite packed matrix using\n*  the Cholesky factorization A = U**T*U or A = L*L**T computed by\n*  SPPTRF.\n*\n*  An estimate is obtained for norm(inv(A)), and the reciprocal of the\n*  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input) REAL array, dimension (N*(N+1)/2)\n*          The triangular factor U or L from the Cholesky factorization\n*          A = U**T*U or A = L*L**T, packed columnwise in a linear\n*          array.  The j-th column of U or L is stored in the array AP\n*          as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.\n*\n*  ANORM   (input) REAL\n*          The 1-norm (or infinity-norm) of the symmetric matrix A.\n*\n*  RCOND   (output) REAL\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an\n*          estimate of the 1-norm of inv(A) computed in this routine.\n*\n*  WORK    (workspace) REAL array, dimension (3*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];
  rb_anorm = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  anorm = (real)NUM2DBL(rb_anorm);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_SFLOAT)
    rb_ap = na_change_type(rb_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rb_ap, real*);
  n = (int)(sqrt((double)8*ldap+1)-1)/2;
  work = ALLOC_N(real, (3*n));
  iwork = ALLOC_N(integer, (n));

  sppcon_(&uplo, &n, ap, &anorm, &rcond, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_sppcon(VALUE mLapack){
  rb_define_module_function(mLapack, "sppcon", rb_sppcon, -1);
}
