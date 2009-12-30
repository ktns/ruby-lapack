#include "rb_lapack.h"

static VALUE
rb_chpevd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_lrwork;
  integer lrwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_w;
  real *w; 
  VALUE rb_z;
  complex *z; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_rwork;
  real *rwork; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  complex *ap_out__;

  integer ldap;
  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, z, work, rwork, iwork, info, ap = NumRu::Lapack.chpevd( jobz, uplo, ap, lwork, lrwork, liwork)\n    or\n  NumRu::Lapack.chpevd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CHPEVD computes all the eigenvalues and, optionally, eigenvectors of\n*  a complex Hermitian matrix A in packed storage.  If eigenvectors are\n*  desired, it uses a divide and conquer algorithm.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) COMPLEX array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*\n*          On exit, AP is overwritten by values generated during the\n*          reduction to tridiagonal form.  If UPLO = 'U', the diagonal\n*          and first superdiagonal of the tridiagonal matrix T overwrite\n*          the corresponding elements of A, and if UPLO = 'L', the\n*          diagonal and first subdiagonal of T overwrite the\n*          corresponding elements of A.\n*\n*  W       (output) REAL array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) COMPLEX array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal\n*          eigenvectors of the matrix A, with the i-th column of Z\n*          holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the required LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of array WORK.\n*          If N <= 1,               LWORK must be at least 1.\n*          If JOBZ = 'N' and N > 1, LWORK must be at least N.\n*          If JOBZ = 'V' and N > 1, LWORK must be at least 2*N.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the required sizes of the WORK, RWORK and\n*          IWORK arrays, returns these values as the first entries of\n*          the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  RWORK   (workspace/output) REAL array, dimension (MAX(1,LRWORK))\n*          On exit, if INFO = 0, RWORK(1) returns the required LRWORK.\n*\n*  LRWORK  (input) INTEGER\n*          The dimension of array RWORK.\n*          If N <= 1,               LRWORK must be at least 1.\n*          If JOBZ = 'N' and N > 1, LRWORK must be at least N.\n*          If JOBZ = 'V' and N > 1, LRWORK must be at least\n*                    1 + 5*N + 2*N**2.\n*\n*          If LRWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the required sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the required LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of array IWORK.\n*          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.\n*          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the required sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, the algorithm failed to converge; i\n*                off-diagonal elements of an intermediate tridiagonal\n*                form did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_jobz = argv[0];
  rb_uplo = argv[1];
  rb_ap = argv[2];
  rb_lwork = argv[3];
  rb_lrwork = argv[4];
  rb_liwork = argv[5];

  jobz = StringValueCStr(rb_jobz)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  lwork = NUM2INT(rb_lwork);
  lrwork = NUM2INT(rb_lrwork);
  liwork = NUM2INT(rb_liwork);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_SCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, complex*);
  n = (int)(sqrt((double)8*ldap+1)-1)/2;
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, complex*);
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
    int shape[1];
    shape[0] = ldap;
    rb_ap_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, complex*);
  MEMCPY(ap_out__, ap, complex, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;

  chpevd_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_w, rb_z, rb_work, rb_rwork, rb_iwork, rb_info, rb_ap);
}

void
init_lapack_chpevd(VALUE mLapack){
  rb_define_module_function(mLapack, "chpevd", rb_chpevd, -1);
}
