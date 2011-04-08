#include "rb_lapack.h"

extern VOID chbgvd_(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, complex *ab, integer *ldab, complex *bb, integer *ldbb, real *w, complex *z, integer *ldz, complex *work, integer *lwork, real *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_chbgvd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobz;
  char jobz; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ka;
  integer ka; 
  VALUE rblapack_kb;
  integer kb; 
  VALUE rblapack_ab;
  complex *ab; 
  VALUE rblapack_bb;
  complex *bb; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_lrwork;
  integer lrwork; 
  VALUE rblapack_liwork;
  integer liwork; 
  VALUE rblapack_w;
  real *w; 
  VALUE rblapack_z;
  complex *z; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_rwork;
  real *rwork; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ab_out__;
  complex *ab_out__;
  VALUE rblapack_bb_out__;
  complex *bb_out__;

  integer ldab;
  integer n;
  integer ldbb;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, z, work, rwork, iwork, info, ab, bb = NumRu::Lapack.chbgvd( jobz, uplo, ka, kb, ab, bb, lwork, lrwork, liwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CHBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CHBGVD computes all the eigenvalues, and optionally, the eigenvectors\n*  of a complex generalized Hermitian-definite banded eigenproblem, of\n*  the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian\n*  and banded, and B is also positive definite.  If eigenvectors are\n*  desired, it uses a divide and conquer algorithm.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangles of A and B are stored;\n*          = 'L':  Lower triangles of A and B are stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  KA      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'. KA >= 0.\n*\n*  KB      (input) INTEGER\n*          The number of superdiagonals of the matrix B if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'. KB >= 0.\n*\n*  AB      (input/output) COMPLEX array, dimension (LDAB, N)\n*          On entry, the upper or lower triangle of the Hermitian band\n*          matrix A, stored in the first ka+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).\n*\n*          On exit, the contents of AB are destroyed.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KA+1.\n*\n*  BB      (input/output) COMPLEX array, dimension (LDBB, N)\n*          On entry, the upper or lower triangle of the Hermitian band\n*          matrix B, stored in the first kb+1 rows of the array.  The\n*          j-th column of B is stored in the j-th column of the array BB\n*          as follows:\n*          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;\n*          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).\n*\n*          On exit, the factor S from the split Cholesky factorization\n*          B = S**H*S, as returned by CPBSTF.\n*\n*  LDBB    (input) INTEGER\n*          The leading dimension of the array BB.  LDBB >= KB+1.\n*\n*  W       (output) REAL array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) COMPLEX array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of\n*          eigenvectors, with the i-th column of Z holding the\n*          eigenvector associated with W(i). The eigenvectors are\n*          normalized so that Z**H*B*Z = I.\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= N.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO=0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If N <= 1,               LWORK >= 1.\n*          If JOBZ = 'N' and N > 1, LWORK >= N.\n*          If JOBZ = 'V' and N > 1, LWORK >= 2*N**2.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal sizes of the WORK, RWORK and\n*          IWORK arrays, returns these values as the first entries of\n*          the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  RWORK   (workspace/output) REAL array, dimension (MAX(1,LRWORK))\n*          On exit, if INFO=0, RWORK(1) returns the optimal LRWORK.\n*\n*  LRWORK  (input) INTEGER\n*          The dimension of array RWORK.\n*          If N <= 1,               LRWORK >= 1.\n*          If JOBZ = 'N' and N > 1, LRWORK >= N.\n*          If JOBZ = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2.\n*\n*          If LRWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO=0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of array IWORK.\n*          If JOBZ = 'N' or N <= 1, LIWORK >= 1.\n*          If JOBZ = 'V' and N > 1, LIWORK >= 3 + 5*N.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is:\n*             <= N:  the algorithm failed to converge:\n*                    i off-diagonal elements of an intermediate\n*                    tridiagonal form did not converge to zero;\n*             > N:   if INFO = N + i, for 1 <= i <= N, then CPBSTF\n*                    returned INFO = i: B is not positive definite.\n*                    The factorization of B could not be completed and\n*                    no eigenvalues or eigenvectors were computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, z, work, rwork, iwork, info, ab, bb = NumRu::Lapack.chbgvd( jobz, uplo, ka, kb, ab, bb, lwork, lrwork, liwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_jobz = argv[0];
  rblapack_uplo = argv[1];
  rblapack_ka = argv[2];
  rblapack_kb = argv[3];
  rblapack_ab = argv[4];
  rblapack_bb = argv[5];
  rblapack_lwork = argv[6];
  rblapack_lrwork = argv[7];
  rblapack_liwork = argv[8];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_bb))
    rb_raise(rb_eArgError, "bb (6th argument) must be NArray");
  if (NA_RANK(rblapack_bb) != 2)
    rb_raise(rb_eArgError, "rank of bb (6th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_bb);
  ldbb = NA_SHAPE0(rblapack_bb);
  if (NA_TYPE(rblapack_bb) != NA_SCOMPLEX)
    rblapack_bb = na_change_type(rblapack_bb, NA_SCOMPLEX);
  bb = NA_PTR_TYPE(rblapack_bb, complex*);
  ka = NUM2INT(rblapack_ka);
  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 1 of bb");
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SCOMPLEX)
    rblapack_ab = na_change_type(rblapack_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rblapack_ab, complex*);
  kb = NUM2INT(rblapack_kb);
  uplo = StringValueCStr(rblapack_uplo)[0];
  jobz = StringValueCStr(rblapack_jobz)[0];
  lrwork = NUM2INT(rblapack_lrwork);
  lwork = NUM2INT(rblapack_lwork);
  liwork = NUM2INT(rblapack_liwork);
  ldz = lsame_(&jobz,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, real*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lrwork);
    rblapack_rwork = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rwork = NA_PTR_TYPE(rblapack_rwork, real*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rblapack_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rblapack_ab_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rblapack_ab_out__, complex*);
  MEMCPY(ab_out__, ab, complex, NA_TOTAL(rblapack_ab));
  rblapack_ab = rblapack_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldbb;
    shape[1] = n;
    rblapack_bb_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  bb_out__ = NA_PTR_TYPE(rblapack_bb_out__, complex*);
  MEMCPY(bb_out__, bb, complex, NA_TOTAL(rblapack_bb));
  rblapack_bb = rblapack_bb_out__;
  bb = bb_out__;

  chbgvd_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(8, rblapack_w, rblapack_z, rblapack_work, rblapack_rwork, rblapack_iwork, rblapack_info, rblapack_ab, rblapack_bb);
}

void
init_lapack_chbgvd(VALUE mLapack){
  rb_define_module_function(mLapack, "chbgvd", rblapack_chbgvd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
