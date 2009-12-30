#include "rb_lapack.h"

static VALUE
rb_cheevd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_lrwork;
  integer lrwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_w;
  real *w; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_rwork;
  real *rwork; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, work, rwork, iwork, info, a = NumRu::Lapack.cheevd( jobz, uplo, a, lwork, lrwork, liwork)\n    or\n  NumRu::Lapack.cheevd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CHEEVD computes all eigenvalues and, optionally, eigenvectors of a\n*  complex Hermitian matrix A.  If eigenvectors are desired, it uses a\n*  divide and conquer algorithm.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA, N)\n*          On entry, the Hermitian matrix A.  If UPLO = 'U', the\n*          leading N-by-N upper triangular part of A contains the\n*          upper triangular part of the matrix A.  If UPLO = 'L',\n*          the leading N-by-N lower triangular part of A contains\n*          the lower triangular part of the matrix A.\n*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the\n*          orthonormal eigenvectors of the matrix A.\n*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')\n*          or the upper triangle (if UPLO='U') of A, including the\n*          diagonal, is destroyed.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  W       (output) REAL array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The length of the array WORK.\n*          If N <= 1,                LWORK must be at least 1.\n*          If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.\n*          If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal sizes of the WORK, RWORK and\n*          IWORK arrays, returns these values as the first entries of\n*          the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  RWORK   (workspace/output) REAL array,\n*                                         dimension (LRWORK)\n*          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.\n*\n*  LRWORK  (input) INTEGER\n*          The dimension of the array RWORK.\n*          If N <= 1,                LRWORK must be at least 1.\n*          If JOBZ  = 'N' and N > 1, LRWORK must be at least N.\n*          If JOBZ  = 'V' and N > 1, LRWORK must be at least\n*                         1 + 5*N + 2*N**2.\n*\n*          If LRWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.\n*          If N <= 1,                LIWORK must be at least 1.\n*          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.\n*          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed\n*                to converge; i off-diagonal elements of an intermediate\n*                tridiagonal form did not converge to zero;\n*                if INFO = i and JOBZ = 'V', then the algorithm failed\n*                to compute an eigenvalue while working on the submatrix\n*                lying in rows and columns INFO/(N+1) through\n*                mod(INFO,N+1).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  Modified description of INFO. Sven, 16 Feb 05.\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_jobz = argv[0];
  rb_uplo = argv[1];
  rb_a = argv[2];
  rb_lwork = argv[3];
  rb_lrwork = argv[4];
  rb_liwork = argv[5];

  jobz = StringValueCStr(rb_jobz)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  lwork = NUM2INT(rb_lwork);
  lrwork = NUM2INT(rb_lrwork);
  liwork = NUM2INT(rb_liwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lrwork);
    rb_rwork = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rwork = NA_PTR_TYPE(rb_rwork, real*);
  {
    int shape[1];
    shape[0] = DIM_LEN(MAX(1,liwork));
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  cheevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_w, rb_work, rb_rwork, rb_iwork, rb_info, rb_a);
}

void
init_lapack_cheevd(VALUE mLapack){
  rb_define_module_function(mLapack, "cheevd", rb_cheevd, -1);
}
