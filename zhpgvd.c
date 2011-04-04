#include "rb_lapack.h"

extern VOID zhpgvd_(integer *itype, char *jobz, char *uplo, integer *n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex *z, integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);

static VALUE
rb_zhpgvd(int argc, VALUE *argv, VALUE self){
  VALUE rb_itype;
  integer itype; 
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_bp;
  doublecomplex *bp; 
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
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  doublecomplex *ap_out__;
  VALUE rb_bp_out__;
  doublecomplex *bp_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer ldap;
  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, z, iwork, info, ap, bp = NumRu::Lapack.zhpgvd( itype, jobz, uplo, ap, bp, lwork, lrwork, liwork)\n    or\n  NumRu::Lapack.zhpgvd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZHPGVD computes all the eigenvalues and, optionally, the eigenvectors\n*  of a complex generalized Hermitian-definite eigenproblem, of the form\n*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and\n*  B are assumed to be Hermitian, stored in packed format, and B is also\n*  positive definite.\n*  If eigenvectors are desired, it uses a divide and conquer algorithm.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  ITYPE   (input) INTEGER\n*          Specifies the problem type to be solved:\n*          = 1:  A*x = (lambda)*B*x\n*          = 2:  A*B*x = (lambda)*x\n*          = 3:  B*A*x = (lambda)*x\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangles of A and B are stored;\n*          = 'L':  Lower triangles of A and B are stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*\n*          On exit, the contents of AP are destroyed.\n*\n*  BP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          B, packed columnwise in a linear array.  The j-th column of B\n*          is stored in the array BP as follows:\n*          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;\n*          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.\n*\n*          On exit, the triangular factor U or L from the Cholesky\n*          factorization B = U**H*U or B = L*L**H, in the same storage\n*          format as B.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of\n*          eigenvectors.  The eigenvectors are normalized as follows:\n*          if ITYPE = 1 or 2, Z**H*B*Z = I;\n*          if ITYPE = 3, Z**H*inv(B)*Z = I.\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the required LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of array WORK.\n*          If N <= 1,               LWORK >= 1.\n*          If JOBZ = 'N' and N > 1, LWORK >= N.\n*          If JOBZ = 'V' and N > 1, LWORK >= 2*N.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the required sizes of the WORK, RWORK and\n*          IWORK arrays, returns these values as the first entries of\n*          the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (MAX(1,LRWORK))\n*          On exit, if INFO = 0, RWORK(1) returns the required LRWORK.\n*\n*  LRWORK  (input) INTEGER\n*          The dimension of array RWORK.\n*          If N <= 1,               LRWORK >= 1.\n*          If JOBZ = 'N' and N > 1, LRWORK >= N.\n*          If JOBZ = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2.\n*\n*          If LRWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the required sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the required LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of array IWORK.\n*          If JOBZ  = 'N' or N <= 1, LIWORK >= 1.\n*          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the required sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  ZPPTRF or ZHPEVD returned an error code:\n*             <= N:  if INFO = i, ZHPEVD failed to converge;\n*                    i off-diagonal elements of an intermediate\n*                    tridiagonal form did not convergeto zero;\n*             > N:   if INFO = N + i, for 1 <= i <= n, then the leading\n*                    minor of order i of B is not positive definite.\n*                    The factorization of B could not be completed and\n*                    no eigenvalues or eigenvectors were computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            LQUERY, UPPER, WANTZ\n      CHARACTER          TRANS\n      INTEGER            J, LIWMIN, LRWMIN, LWMIN, NEIG\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA, ZHPEVD, ZHPGST, ZPPTRF, ZTPMV, ZTPSV\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          DBLE, MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_itype = argv[0];
  rb_jobz = argv[1];
  rb_uplo = argv[2];
  rb_ap = argv[3];
  rb_bp = argv[4];
  rb_lwork = argv[5];
  rb_lrwork = argv[6];
  rb_liwork = argv[7];

  uplo = StringValueCStr(rb_uplo)[0];
  jobz = StringValueCStr(rb_jobz)[0];
  liwork = NUM2INT(rb_liwork);
  lrwork = NUM2INT(rb_lrwork);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  itype = NUM2INT(rb_itype);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  if (!NA_IsNArray(rb_bp))
    rb_raise(rb_eArgError, "bp (5th argument) must be NArray");
  if (NA_RANK(rb_bp) != 1)
    rb_raise(rb_eArgError, "rank of bp (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_bp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of bp must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_bp) != NA_DCOMPLEX)
    rb_bp = na_change_type(rb_bp, NA_DCOMPLEX);
  bp = NA_PTR_TYPE(rb_bp, doublecomplex*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  {
    int shape[1];
    shape[0] = ldap;
    rb_ap_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, doublecomplex*);
  MEMCPY(ap_out__, ap, doublecomplex, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_bp_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  bp_out__ = NA_PTR_TYPE(rb_bp_out__, doublecomplex*);
  MEMCPY(bp_out__, bp, doublecomplex, NA_TOTAL(rb_bp));
  rb_bp = rb_bp_out__;
  bp = bp_out__;
  work = ALLOC_N(doublecomplex, (MAX(1,lwork)));
  rwork = ALLOC_N(doublereal, (MAX(1,lrwork)));

  zhpgvd_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  free(work);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_w, rb_z, rb_iwork, rb_info, rb_ap, rb_bp);
}

void
init_lapack_zhpgvd(VALUE mLapack){
  rb_define_module_function(mLapack, "zhpgvd", rb_zhpgvd, -1);
}
