#include "rb_lapack.h"

static VALUE
rb_zhbgvd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ka;
  integer ka; 
  VALUE rb_kb;
  integer kb; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb_bb;
  doublecomplex *bb; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_lrwork;
  integer lrwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublecomplex *z; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_rwork;
  doublereal *rwork; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublecomplex *ab_out__;
  VALUE rb_bb_out__;
  doublecomplex *bb_out__;

  integer ldab;
  integer n;
  integer ldbb;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, z, work, rwork, iwork, info, ab, bb = NumRu::Lapack.zhbgvd( jobz, uplo, ka, kb, ab, bb, lwork, lrwork, liwork)\n    or\n  NumRu::Lapack.zhbgvd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZHBGVD computes all the eigenvalues, and optionally, the eigenvectors\n*  of a complex generalized Hermitian-definite banded eigenproblem, of\n*  the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian\n*  and banded, and B is also positive definite.  If eigenvectors are\n*  desired, it uses a divide and conquer algorithm.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangles of A and B are stored;\n*          = 'L':  Lower triangles of A and B are stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  KA      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'. KA >= 0.\n*\n*  KB      (input) INTEGER\n*          The number of superdiagonals of the matrix B if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'. KB >= 0.\n*\n*  AB      (input/output) COMPLEX*16 array, dimension (LDAB, N)\n*          On entry, the upper or lower triangle of the Hermitian band\n*          matrix A, stored in the first ka+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).\n*\n*          On exit, the contents of AB are destroyed.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KA+1.\n*\n*  BB      (input/output) COMPLEX*16 array, dimension (LDBB, N)\n*          On entry, the upper or lower triangle of the Hermitian band\n*          matrix B, stored in the first kb+1 rows of the array.  The\n*          j-th column of B is stored in the j-th column of the array BB\n*          as follows:\n*          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;\n*          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).\n*\n*          On exit, the factor S from the split Cholesky factorization\n*          B = S**H*S, as returned by ZPBSTF.\n*\n*  LDBB    (input) INTEGER\n*          The leading dimension of the array BB.  LDBB >= KB+1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of\n*          eigenvectors, with the i-th column of Z holding the\n*          eigenvector associated with W(i). The eigenvectors are\n*          normalized so that Z**H*B*Z = I.\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= N.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO=0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If N <= 1,               LWORK >= 1.\n*          If JOBZ = 'N' and N > 1, LWORK >= N.\n*          If JOBZ = 'V' and N > 1, LWORK >= 2*N**2.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal sizes of the WORK, RWORK and\n*          IWORK arrays, returns these values as the first entries of\n*          the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  RWORK   (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LRWORK))\n*          On exit, if INFO=0, RWORK(1) returns the optimal LRWORK.\n*\n*  LRWORK  (input) INTEGER\n*          The dimension of array RWORK.\n*          If N <= 1,               LRWORK >= 1.\n*          If JOBZ = 'N' and N > 1, LRWORK >= N.\n*          If JOBZ = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2.\n*\n*          If LRWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO=0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of array IWORK.\n*          If JOBZ = 'N' or N <= 1, LIWORK >= 1.\n*          If JOBZ = 'V' and N > 1, LIWORK >= 3 + 5*N.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is:\n*             <= N:  the algorithm failed to converge:\n*                    i off-diagonal elements of an intermediate\n*                    tridiagonal form did not converge to zero;\n*             > N:   if INFO = N + i, for 1 <= i <= N, then ZPBSTF\n*                    returned INFO = i: B is not positive definite.\n*                    The factorization of B could not be completed and\n*                    no eigenvalues or eigenvectors were computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_jobz = argv[0];
  rb_uplo = argv[1];
  rb_ka = argv[2];
  rb_kb = argv[3];
  rb_ab = argv[4];
  rb_bb = argv[5];
  rb_lwork = argv[6];
  rb_lrwork = argv[7];
  rb_liwork = argv[8];

  jobz = StringValueCStr(rb_jobz)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  ka = NUM2INT(rb_ka);
  kb = NUM2INT(rb_kb);
  lwork = NUM2INT(rb_lwork);
  lrwork = NUM2INT(rb_lrwork);
  liwork = NUM2INT(rb_liwork);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  if (!NA_IsNArray(rb_bb))
    rb_raise(rb_eArgError, "bb (6th argument) must be NArray");
  if (NA_RANK(rb_bb) != 2)
    rb_raise(rb_eArgError, "rank of bb (6th argument) must be %d", 2);
  ldbb = NA_SHAPE0(rb_bb);
  if (NA_SHAPE1(rb_bb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of bb must be the same as shape 1 of ab");
  if (NA_TYPE(rb_bb) != NA_DCOMPLEX)
    rb_bb = na_change_type(rb_bb, NA_DCOMPLEX);
  bb = NA_PTR_TYPE(rb_bb, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  ldz = lsame_(&jobz,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lrwork);
    rb_rwork = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  rwork = NA_PTR_TYPE(rb_rwork, doublereal*);
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
    rb_ab_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublecomplex*);
  MEMCPY(ab_out__, ab, doublecomplex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldbb;
    shape[1] = n;
    rb_bb_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  bb_out__ = NA_PTR_TYPE(rb_bb_out__, doublecomplex*);
  MEMCPY(bb_out__, bb, doublecomplex, NA_TOTAL(rb_bb));
  rb_bb = rb_bb_out__;
  bb = bb_out__;

  zhbgvd_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_w, rb_z, rb_work, rb_rwork, rb_iwork, rb_info, rb_ab, rb_bb);
}

void
init_lapack_zhbgvd(VALUE mLapack){
  rb_define_module_function(mLapack, "zhbgvd", rb_zhbgvd, -1);
}
