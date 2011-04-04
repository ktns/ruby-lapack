#include "rb_lapack.h"

extern VOID ssbevd_(char *jobz, char *uplo, integer *n, integer *kd, real *ab, integer *ldab, real *w, real *z, integer *ldz, real *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE
rb_ssbevd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_w;
  real *w; 
  VALUE rb_z;
  real *z; 
  VALUE rb_work;
  real *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  real *ab_out__;

  integer ldab;
  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, z, work, iwork, info, ab = NumRu::Lapack.ssbevd( jobz, uplo, kd, ab, lwork, liwork)\n    or\n  NumRu::Lapack.ssbevd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSBEVD computes all the eigenvalues and, optionally, eigenvectors of\n*  a real symmetric band matrix A. If eigenvectors are desired, it uses\n*  a divide and conquer algorithm.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input/output) REAL array, dimension (LDAB, N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first KD+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*\n*          On exit, AB is overwritten by values generated during the\n*          reduction to tridiagonal form.  If UPLO = 'U', the first\n*          superdiagonal and the diagonal of the tridiagonal matrix T\n*          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',\n*          the diagonal and first subdiagonal of T are returned in the\n*          first two rows of AB.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD + 1.\n*\n*  W       (output) REAL array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) REAL array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal\n*          eigenvectors of the matrix A, with the i-th column of Z\n*          holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace/output) REAL array,\n*                                         dimension (LWORK)\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          IF N <= 1,                LWORK must be at least 1.\n*          If JOBZ  = 'N' and N > 2, LWORK must be at least 2*N.\n*          If JOBZ  = 'V' and N > 2, LWORK must be at least\n*                         ( 1 + 5*N + 2*N**2 ).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal sizes of the WORK and IWORK\n*          arrays, returns these values as the first entries of the WORK\n*          and IWORK arrays, and no error message related to LWORK or\n*          LIWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array LIWORK.\n*          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.\n*          If JOBZ  = 'V' and N > 2, LIWORK must be at least 3 + 5*N.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK and\n*          IWORK arrays, returns these values as the first entries of\n*          the WORK and IWORK arrays, and no error message related to\n*          LWORK or LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the algorithm failed to converge; i\n*                off-diagonal elements of an intermediate tridiagonal\n*                form did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_jobz = argv[0];
  rb_uplo = argv[1];
  rb_kd = argv[2];
  rb_ab = argv[3];
  rb_lwork = argv[4];
  rb_liwork = argv[5];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  liwork = NUM2INT(rb_liwork);
  jobz = StringValueCStr(rb_jobz)[0];
  lwork = NUM2INT(rb_lwork);
  kd = NUM2INT(rb_kd);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;

  ssbevd_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, iwork, &liwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_w, rb_z, rb_work, rb_iwork, rb_info, rb_ab);
}

void
init_lapack_ssbevd(VALUE mLapack){
  rb_define_module_function(mLapack, "ssbevd", rb_ssbevd, -1);
}
