#include "rb_lapack.h"

static VALUE
rb_zhpgvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_itype;
  integer itype; 
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_bp;
  doublecomplex *bp; 
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_abstol;
  doublereal abstol; 
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublecomplex *z; 
  VALUE rb_ifail;
  integer *ifail; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  doublecomplex *ap_out__;
  VALUE rb_bp_out__;
  doublecomplex *bp_out__;
  doublecomplex *work;
  doublereal *rwork;
  integer *iwork;

  integer ldap;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, w, z, ifail, info, ap, bp = NumRu::Lapack.zhpgvx( itype, jobz, range, uplo, ap, bp, vl, vu, il, iu, abstol, ldz)\n    or\n  NumRu::Lapack.zhpgvx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZHPGVX computes selected eigenvalues and, optionally, eigenvectors\n*  of a complex generalized Hermitian-definite eigenproblem, of the form\n*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and\n*  B are assumed to be Hermitian, stored in packed format, and B is also\n*  positive definite.  Eigenvalues and eigenvectors can be selected by\n*  specifying either a range of values or a range of indices for the\n*  desired eigenvalues.\n*\n\n*  Arguments\n*  =========\n*\n*  ITYPE   (input) INTEGER\n*          Specifies the problem type to be solved:\n*          = 1:  A*x = (lambda)*B*x\n*          = 2:  A*B*x = (lambda)*x\n*          = 3:  B*A*x = (lambda)*x\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found;\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found;\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangles of A and B are stored;\n*          = 'L':  Lower triangles of A and B are stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*\n*          On exit, the contents of AP are destroyed.\n*\n*  BP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          B, packed columnwise in a linear array.  The j-th column of B\n*          is stored in the array BP as follows:\n*          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;\n*          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.\n*\n*          On exit, the triangular factor U or L from the Cholesky\n*          factorization B = U**H*U or B = L*L**H, in the same storage\n*          format as B.\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) DOUBLE PRECISION\n*          The absolute error tolerance for the eigenvalues.\n*          An approximate eigenvalue is accepted as converged\n*          when it is determined to lie in an interval [a,b]\n*          of width less than or equal to\n*\n*                  ABSTOL + EPS *   max( |a|,|b| ) ,\n*\n*          where EPS is the machine precision.  If ABSTOL is less than\n*          or equal to zero, then  EPS*|T|  will be used in its place,\n*          where |T| is the 1-norm of the tridiagonal matrix obtained\n*          by reducing AP to tridiagonal form.\n*\n*          Eigenvalues will be computed most accurately when ABSTOL is\n*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.\n*          If this routine returns with INFO>0, indicating that some\n*          eigenvectors did not converge, try setting ABSTOL to\n*          2*DLAMCH('S').\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          On normal exit, the first M elements contain the selected\n*          eigenvalues in ascending order.\n*\n*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)\n*          If JOBZ = 'N', then Z is not referenced.\n*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix A\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          The eigenvectors are normalized as follows:\n*          if ITYPE = 1 or 2, Z**H*B*Z = I;\n*          if ITYPE = 3, Z**H*inv(B)*Z = I.\n*\n*          If an eigenvector fails to converge, then that column of Z\n*          contains the latest approximation to the eigenvector, and the\n*          index of the eigenvector is returned in IFAIL.\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and an upper bound must be used.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (2*N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (7*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (5*N)\n*\n*  IFAIL   (output) INTEGER array, dimension (N)\n*          If JOBZ = 'V', then if INFO = 0, the first M elements of\n*          IFAIL are zero.  If INFO > 0, then IFAIL contains the\n*          indices of the eigenvectors that failed to converge.\n*          If JOBZ = 'N', then IFAIL is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  ZPPTRF or ZHPEVX returned an error code:\n*             <= N:  if INFO = i, ZHPEVX failed to converge;\n*                    i eigenvectors failed to converge.  Their indices\n*                    are stored in array IFAIL.\n*             > N:   if INFO = N + i, for 1 <= i <= n, then the leading\n*                    minor of order i of B is not positive definite.\n*                    The factorization of B could not be completed and\n*                    no eigenvalues or eigenvectors were computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            ALLEIG, INDEIG, UPPER, VALEIG, WANTZ\n      CHARACTER          TRANS\n      INTEGER            J\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA, ZHPEVX, ZHPGST, ZPPTRF, ZTPMV, ZTPSV\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_itype = argv[0];
  rb_jobz = argv[1];
  rb_range = argv[2];
  rb_uplo = argv[3];
  rb_ap = argv[4];
  rb_bp = argv[5];
  rb_vl = argv[6];
  rb_vu = argv[7];
  rb_il = argv[8];
  rb_iu = argv[9];
  rb_abstol = argv[10];
  rb_ldz = argv[11];

  itype = NUM2INT(rb_itype);
  jobz = StringValueCStr(rb_jobz)[0];
  range = StringValueCStr(rb_range)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  vl = NUM2DBL(rb_vl);
  vu = NUM2DBL(rb_vu);
  il = NUM2INT(rb_il);
  iu = NUM2INT(rb_iu);
  abstol = NUM2DBL(rb_abstol);
  ldz = NUM2INT(rb_ldz);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  if (!NA_IsNArray(rb_bp))
    rb_raise(rb_eArgError, "bp (6th argument) must be NArray");
  if (NA_RANK(rb_bp) != 1)
    rb_raise(rb_eArgError, "rank of bp (6th argument) must be %d", 1);
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
    shape[0] = lsame_(&jobz,"N") ? 0 : ldz;
    shape[1] = lsame_(&jobz,"N") ? 0 : n;
    rb_z = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_ifail = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifail = NA_PTR_TYPE(rb_ifail, integer*);
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
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (7*n));
  iwork = ALLOC_N(integer, (5*n));

  zhpgvx_(&itype, &jobz, &range, &uplo, &n, ap, bp, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, rwork, iwork, ifail, &info);

  free(work);
  free(rwork);
  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_m, rb_w, rb_z, rb_ifail, rb_info, rb_ap, rb_bp);
}

void
init_lapack_zhpgvx(VALUE mLapack){
  rb_define_module_function(mLapack, "zhpgvx", rb_zhpgvx, -1);
}
