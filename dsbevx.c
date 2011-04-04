#include "rb_lapack.h"

extern VOID dsbevx_(char *jobz, char *range, char *uplo, integer *n, integer *kd, doublereal *ab, integer *ldab, doublereal *q, integer *ldq, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z, integer *ldz, doublereal *work, integer *iwork, integer *ifail, integer *info);

static VALUE
rb_dsbevx(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  doublereal *ab; 
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
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_ifail;
  integer *ifail; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublereal *ab_out__;
  doublereal *work;
  integer *iwork;

  integer ldab;
  integer n;
  integer ldq;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  q, m, w, z, ifail, info, ab = NumRu::Lapack.dsbevx( jobz, range, uplo, kd, ab, vl, vu, il, iu, abstol)\n    or\n  NumRu::Lapack.dsbevx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSBEVX computes selected eigenvalues and, optionally, eigenvectors\n*  of a real symmetric band matrix A.  Eigenvalues and eigenvectors can\n*  be selected by specifying either a range of values or a range of\n*  indices for the desired eigenvalues.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found;\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found;\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first KD+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*\n*          On exit, AB is overwritten by values generated during the\n*          reduction to tridiagonal form.  If UPLO = 'U', the first\n*          superdiagonal and the diagonal of the tridiagonal matrix T\n*          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',\n*          the diagonal and first subdiagonal of T are returned in the\n*          first two rows of AB.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD + 1.\n*\n*  Q       (output) DOUBLE PRECISION array, dimension (LDQ, N)\n*          If JOBZ = 'V', the N-by-N orthogonal matrix used in the\n*                         reduction to tridiagonal form.\n*          If JOBZ = 'N', the array Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  If JOBZ = 'V', then\n*          LDQ >= max(1,N).\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) DOUBLE PRECISION\n*          The absolute error tolerance for the eigenvalues.\n*          An approximate eigenvalue is accepted as converged\n*          when it is determined to lie in an interval [a,b]\n*          of width less than or equal to\n*\n*                  ABSTOL + EPS *   max( |a|,|b| ) ,\n*\n*          where EPS is the machine precision.  If ABSTOL is less than\n*          or equal to zero, then  EPS*|T|  will be used in its place,\n*          where |T| is the 1-norm of the tridiagonal matrix obtained\n*          by reducing AB to tridiagonal form.\n*\n*          Eigenvalues will be computed most accurately when ABSTOL is\n*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.\n*          If this routine returns with INFO>0, indicating that some\n*          eigenvectors did not converge, try setting ABSTOL to\n*          2*DLAMCH('S').\n*\n*          See \"Computing Small Singular Values of Bidiagonal Matrices\n*          with Guaranteed High Relative Accuracy,\" by Demmel and\n*          Kahan, LAPACK Working Note #3.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          The first M elements contain the selected eigenvalues in\n*          ascending order.\n*\n*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))\n*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix A\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          If an eigenvector fails to converge, then that column of Z\n*          contains the latest approximation to the eigenvector, and the\n*          index of the eigenvector is returned in IFAIL.\n*          If JOBZ = 'N', then Z is not referenced.\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and an upper bound must be used.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (7*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (5*N)\n*\n*  IFAIL   (output) INTEGER array, dimension (N)\n*          If JOBZ = 'V', then if INFO = 0, the first M elements of\n*          IFAIL are zero.  If INFO > 0, then IFAIL contains the\n*          indices of the eigenvectors that failed to converge.\n*          If JOBZ = 'N', then IFAIL is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, then i eigenvectors failed to converge.\n*                Their indices are stored in array IFAIL.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_jobz = argv[0];
  rb_range = argv[1];
  rb_uplo = argv[2];
  rb_kd = argv[3];
  rb_ab = argv[4];
  rb_vl = argv[5];
  rb_vu = argv[6];
  rb_il = argv[7];
  rb_iu = argv[8];
  rb_abstol = argv[9];

  abstol = NUM2DBL(rb_abstol);
  vl = NUM2DBL(rb_vl);
  iu = NUM2INT(rb_iu);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DFLOAT)
    rb_ab = na_change_type(rb_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rb_ab, doublereal*);
  jobz = StringValueCStr(rb_jobz)[0];
  vu = NUM2DBL(rb_vu);
  range = StringValueCStr(rb_range)[0];
  il = NUM2INT(rb_il);
  kd = NUM2INT(rb_kd);
  uplo = StringValueCStr(rb_uplo)[0];
  m = lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0;
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  ldq = lsame_(&jobz,"V") ? MAX(1,n) : 0;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rb_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_ifail = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifail = NA_PTR_TYPE(rb_ifail, integer*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublereal*);
  MEMCPY(ab_out__, ab, doublereal, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  work = ALLOC_N(doublereal, (7*n));
  iwork = ALLOC_N(integer, (5*n));

  dsbevx_(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, iwork, ifail, &info);

  free(work);
  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_q, rb_m, rb_w, rb_z, rb_ifail, rb_info, rb_ab);
}

void
init_lapack_dsbevx(VALUE mLapack){
  rb_define_module_function(mLapack, "dsbevx", rb_dsbevx, -1);
}
