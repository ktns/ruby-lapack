#include "rb_lapack.h"

static VALUE
rb_zhpgv(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublecomplex *z; 
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
    printf("%s\n", "USAGE:\n  w, z, info, ap, bp = NumRu::Lapack.zhpgv( itype, jobz, uplo, ap, bp)\n    or\n  NumRu::Lapack.zhpgv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZHPGV computes all the eigenvalues and, optionally, the eigenvectors\n*  of a complex generalized Hermitian-definite eigenproblem, of the form\n*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.\n*  Here A and B are assumed to be Hermitian, stored in packed format,\n*  and B is also positive definite.\n*\n\n*  Arguments\n*  =========\n*\n*  ITYPE   (input) INTEGER\n*          Specifies the problem type to be solved:\n*          = 1:  A*x = (lambda)*B*x\n*          = 2:  A*B*x = (lambda)*x\n*          = 3:  B*A*x = (lambda)*x\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangles of A and B are stored;\n*          = 'L':  Lower triangles of A and B are stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*\n*          On exit, the contents of AP are destroyed.\n*\n*  BP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          B, packed columnwise in a linear array.  The j-th column of B\n*          is stored in the array BP as follows:\n*          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;\n*          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.\n*\n*          On exit, the triangular factor U or L from the Cholesky\n*          factorization B = U**H*U or B = L*L**H, in the same storage\n*          format as B.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of\n*          eigenvectors.  The eigenvectors are normalized as follows:\n*          if ITYPE = 1 or 2, Z**H*B*Z = I;\n*          if ITYPE = 3, Z**H*inv(B)*Z = I.\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (max(1, 2*N-1))\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2))\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  ZPPTRF or ZHPEV returned an error code:\n*             <= N:  if INFO = i, ZHPEV failed to converge;\n*                    i off-diagonal elements of an intermediate\n*                    tridiagonal form did not convergeto zero;\n*             > N:   if INFO = N + i, for 1 <= i <= n, then the leading\n*                    minor of order i of B is not positive definite.\n*                    The factorization of B could not be completed and\n*                    no eigenvalues or eigenvectors were computed.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            UPPER, WANTZ\n      CHARACTER          TRANS\n      INTEGER            J, NEIG\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA, ZHPEV, ZHPGST, ZPPTRF, ZTPMV, ZTPSV\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_itype = argv[0];
  rb_jobz = argv[1];
  rb_uplo = argv[2];
  rb_ap = argv[3];
  rb_bp = argv[4];

  itype = NUM2INT(rb_itype);
  jobz = StringValueCStr(rb_jobz)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  if (!NA_IsNArray(rb_bp))
    rb_raise(rb_eArgError, "bp (5th argument) must be NArray");
  if (NA_RANK(rb_bp) != 1)
    rb_raise(rb_eArgError, "rank of bp (5th argument) must be %d", 1);
  n = (int)(sqrt((double)8*ldap+1)-1)/2;
  if (NA_SHAPE0(rb_bp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of bp must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_bp) != NA_DCOMPLEX)
    rb_bp = na_change_type(rb_bp, NA_DCOMPLEX);
  bp = NA_PTR_TYPE(rb_bp, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
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
  work = ALLOC_N(doublecomplex, (MAX(1, 2*n-1)));
  rwork = ALLOC_N(doublereal, (MAX(1, 3*n-2)));

  zhpgv_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, rwork, &info);

  free(work);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_w, rb_z, rb_info, rb_ap, rb_bp);
}

void
init_lapack_zhpgv(VALUE mLapack){
  rb_define_module_function(mLapack, "zhpgv", rb_zhpgv, -1);
}
