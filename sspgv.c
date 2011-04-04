#include "rb_lapack.h"

extern VOID sspgv_(integer *itype, char *jobz, char *uplo, integer *n, real *ap, real *bp, real *w, real *z, integer *ldz, real *work, integer *info);

static VALUE
rb_sspgv(int argc, VALUE *argv, VALUE self){
  VALUE rb_itype;
  integer itype; 
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_bp;
  real *bp; 
  VALUE rb_w;
  real *w; 
  VALUE rb_z;
  real *z; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  real *ap_out__;
  VALUE rb_bp_out__;
  real *bp_out__;
  real *work;

  integer ldap;
  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, z, info, ap, bp = NumRu::Lapack.sspgv( itype, jobz, uplo, ap, bp)\n    or\n  NumRu::Lapack.sspgv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSPGV computes all the eigenvalues and, optionally, the eigenvectors\n*  of a real generalized symmetric-definite eigenproblem, of the form\n*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.\n*  Here A and B are assumed to be symmetric, stored in packed format,\n*  and B is also positive definite.\n*\n\n*  Arguments\n*  =========\n*\n*  ITYPE   (input) INTEGER\n*          Specifies the problem type to be solved:\n*          = 1:  A*x = (lambda)*B*x\n*          = 2:  A*B*x = (lambda)*x\n*          = 3:  B*A*x = (lambda)*x\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangles of A and B are stored;\n*          = 'L':  Lower triangles of A and B are stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  AP      (input/output) REAL array, dimension\n*                            (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the symmetric matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*\n*          On exit, the contents of AP are destroyed.\n*\n*  BP      (input/output) REAL array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the symmetric matrix\n*          B, packed columnwise in a linear array.  The j-th column of B\n*          is stored in the array BP as follows:\n*          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;\n*          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.\n*\n*          On exit, the triangular factor U or L from the Cholesky\n*          factorization B = U**T*U or B = L*L**T, in the same storage\n*          format as B.\n*\n*  W       (output) REAL array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) REAL array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of\n*          eigenvectors.  The eigenvectors are normalized as follows:\n*          if ITYPE = 1 or 2, Z**T*B*Z = I;\n*          if ITYPE = 3, Z**T*inv(B)*Z = I.\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (3*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  SPPTRF or SSPEV returned an error code:\n*             <= N:  if INFO = i, SSPEV failed to converge;\n*                    i off-diagonal elements of an intermediate\n*                    tridiagonal form did not converge to zero.\n*             > N:   if INFO = n + i, for 1 <= i <= n, then the leading\n*                    minor of order i of B is not positive definite.\n*                    The factorization of B could not be completed and\n*                    no eigenvalues or eigenvectors were computed.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            UPPER, WANTZ\n      CHARACTER          TRANS\n      INTEGER            J, NEIG\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SPPTRF, SSPEV, SSPGST, STPMV, STPSV, XERBLA\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_itype = argv[0];
  rb_jobz = argv[1];
  rb_uplo = argv[2];
  rb_ap = argv[3];
  rb_bp = argv[4];

  uplo = StringValueCStr(rb_uplo)[0];
  jobz = StringValueCStr(rb_jobz)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_SFLOAT)
    rb_ap = na_change_type(rb_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rb_ap, real*);
  itype = NUM2INT(rb_itype);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  if (!NA_IsNArray(rb_bp))
    rb_raise(rb_eArgError, "bp (5th argument) must be NArray");
  if (NA_RANK(rb_bp) != 1)
    rb_raise(rb_eArgError, "rank of bp (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_bp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of bp must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_bp) != NA_SFLOAT)
    rb_bp = na_change_type(rb_bp, NA_SFLOAT);
  bp = NA_PTR_TYPE(rb_bp, real*);
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
    shape[0] = ldap;
    rb_ap_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, real*);
  MEMCPY(ap_out__, ap, real, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_bp_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  bp_out__ = NA_PTR_TYPE(rb_bp_out__, real*);
  MEMCPY(bp_out__, bp, real, NA_TOTAL(rb_bp));
  rb_bp = rb_bp_out__;
  bp = bp_out__;
  work = ALLOC_N(real, (3*n));

  sspgv_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_w, rb_z, rb_info, rb_ap, rb_bp);
}

void
init_lapack_sspgv(VALUE mLapack){
  rb_define_module_function(mLapack, "sspgv", rb_sspgv, -1);
}
