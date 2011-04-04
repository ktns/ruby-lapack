#include "rb_lapack.h"

extern VOID dpbcon_(char *uplo, integer *n, integer *kd, doublereal *ab, integer *ldab, doublereal *anorm, doublereal *rcond, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dpbcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  doublereal *ab; 
  VALUE rb_anorm;
  doublereal anorm; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_info;
  integer info; 
  doublereal *work;
  integer *iwork;

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.dpbcon( uplo, kd, ab, anorm)\n    or\n  NumRu::Lapack.dpbcon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DPBCON estimates the reciprocal of the condition number (in the\n*  1-norm) of a real symmetric positive definite band matrix using the\n*  Cholesky factorization A = U**T*U or A = L*L**T computed by DPBTRF.\n*\n*  An estimate is obtained for norm(inv(A)), and the reciprocal of the\n*  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangular factor stored in AB;\n*          = 'L':  Lower triangular factor stored in AB.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)\n*          The triangular factor U or L from the Cholesky factorization\n*          A = U**T*U or A = L*L**T of the band matrix A, stored in the\n*          first KD+1 rows of the array.  The j-th column of U or L is\n*          stored in the j-th column of the array AB as follows:\n*          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD+1.\n*\n*  ANORM   (input) DOUBLE PRECISION\n*          The 1-norm (or infinity-norm) of the symmetric band matrix A.\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an\n*          estimate of the 1-norm of inv(A) computed in this routine.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_kd = argv[1];
  rb_ab = argv[2];
  rb_anorm = argv[3];

  anorm = NUM2DBL(rb_anorm);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DFLOAT)
    rb_ab = na_change_type(rb_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rb_ab, doublereal*);
  kd = NUM2INT(rb_kd);
  uplo = StringValueCStr(rb_uplo)[0];
  work = ALLOC_N(doublereal, (3*n));
  iwork = ALLOC_N(integer, (n));

  dpbcon_(&uplo, &n, &kd, ab, &ldab, &anorm, &rcond, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_dpbcon(VALUE mLapack){
  rb_define_module_function(mLapack, "dpbcon", rb_dpbcon, -1);
}
