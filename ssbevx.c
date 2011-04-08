#include "rb_lapack.h"

extern VOID ssbevx_(char *jobz, char *range, char *uplo, integer *n, integer *kd, real *ab, integer *ldab, real *q, integer *ldq, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, real *z, integer *ldz, real *work, integer *iwork, integer *ifail, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ssbevx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobz;
  char jobz; 
  VALUE rblapack_range;
  char range; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_kd;
  integer kd; 
  VALUE rblapack_ab;
  real *ab; 
  VALUE rblapack_vl;
  real vl; 
  VALUE rblapack_vu;
  real vu; 
  VALUE rblapack_il;
  integer il; 
  VALUE rblapack_iu;
  integer iu; 
  VALUE rblapack_abstol;
  real abstol; 
  VALUE rblapack_q;
  real *q; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_w;
  real *w; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_ifail;
  integer *ifail; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ab_out__;
  real *ab_out__;
  real *work;
  integer *iwork;

  integer ldab;
  integer n;
  integer ldq;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  q, m, w, z, ifail, info, ab = NumRu::Lapack.ssbevx( jobz, range, uplo, kd, ab, vl, vu, il, iu, abstol, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSBEVX computes selected eigenvalues and, optionally, eigenvectors\n*  of a real symmetric band matrix A.  Eigenvalues and eigenvectors can\n*  be selected by specifying either a range of values or a range of\n*  indices for the desired eigenvalues.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found;\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found;\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input/output) REAL array, dimension (LDAB, N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first KD+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*\n*          On exit, AB is overwritten by values generated during the\n*          reduction to tridiagonal form.  If UPLO = 'U', the first\n*          superdiagonal and the diagonal of the tridiagonal matrix T\n*          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',\n*          the diagonal and first subdiagonal of T are returned in the\n*          first two rows of AB.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD + 1.\n*\n*  Q       (output) REAL array, dimension (LDQ, N)\n*          If JOBZ = 'V', the N-by-N orthogonal matrix used in the\n*                         reduction to tridiagonal form.\n*          If JOBZ = 'N', the array Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  If JOBZ = 'V', then\n*          LDQ >= max(1,N).\n*\n*  VL      (input) REAL\n*  VU      (input) REAL\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) REAL\n*          The absolute error tolerance for the eigenvalues.\n*          An approximate eigenvalue is accepted as converged\n*          when it is determined to lie in an interval [a,b]\n*          of width less than or equal to\n*\n*                  ABSTOL + EPS *   max( |a|,|b| ) ,\n*\n*          where EPS is the machine precision.  If ABSTOL is less than\n*          or equal to zero, then  EPS*|T|  will be used in its place,\n*          where |T| is the 1-norm of the tridiagonal matrix obtained\n*          by reducing AB to tridiagonal form.\n*\n*          Eigenvalues will be computed most accurately when ABSTOL is\n*          set to twice the underflow threshold 2*SLAMCH('S'), not zero.\n*          If this routine returns with INFO>0, indicating that some\n*          eigenvectors did not converge, try setting ABSTOL to\n*          2*SLAMCH('S').\n*\n*          See \"Computing Small Singular Values of Bidiagonal Matrices\n*          with Guaranteed High Relative Accuracy,\" by Demmel and\n*          Kahan, LAPACK Working Note #3.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) REAL array, dimension (N)\n*          The first M elements contain the selected eigenvalues in\n*          ascending order.\n*\n*  Z       (output) REAL array, dimension (LDZ, max(1,M))\n*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix A\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          If an eigenvector fails to converge, then that column of Z\n*          contains the latest approximation to the eigenvector, and the\n*          index of the eigenvector is returned in IFAIL.\n*          If JOBZ = 'N', then Z is not referenced.\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and an upper bound must be used.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (7*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (5*N)\n*\n*  IFAIL   (output) INTEGER array, dimension (N)\n*          If JOBZ = 'V', then if INFO = 0, the first M elements of\n*          IFAIL are zero.  If INFO > 0, then IFAIL contains the\n*          indices of the eigenvectors that failed to converge.\n*          If JOBZ = 'N', then IFAIL is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, then i eigenvectors failed to converge.\n*                Their indices are stored in array IFAIL.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  q, m, w, z, ifail, info, ab = NumRu::Lapack.ssbevx( jobz, range, uplo, kd, ab, vl, vu, il, iu, abstol, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_jobz = argv[0];
  rblapack_range = argv[1];
  rblapack_uplo = argv[2];
  rblapack_kd = argv[3];
  rblapack_ab = argv[4];
  rblapack_vl = argv[5];
  rblapack_vu = argv[6];
  rblapack_il = argv[7];
  rblapack_iu = argv[8];
  rblapack_abstol = argv[9];
  if (rb_options != Qnil) {
  }

  abstol = (real)NUM2DBL(rblapack_abstol);
  vl = (real)NUM2DBL(rblapack_vl);
  iu = NUM2INT(rblapack_iu);
  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_ab);
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SFLOAT)
    rblapack_ab = na_change_type(rblapack_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rblapack_ab, real*);
  jobz = StringValueCStr(rblapack_jobz)[0];
  vu = (real)NUM2DBL(rblapack_vu);
  range = StringValueCStr(rblapack_range)[0];
  il = NUM2INT(rblapack_il);
  kd = NUM2INT(rblapack_kd);
  uplo = StringValueCStr(rblapack_uplo)[0];
  m = lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0;
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  ldq = lsame_(&jobz,"V") ? MAX(1,n) : 0;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rblapack_q, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, real*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rblapack_z = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_ifail = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifail = NA_PTR_TYPE(rblapack_ifail, integer*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rblapack_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rblapack_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rblapack_ab));
  rblapack_ab = rblapack_ab_out__;
  ab = ab_out__;
  work = ALLOC_N(real, (7*n));
  iwork = ALLOC_N(integer, (5*n));

  ssbevx_(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, iwork, ifail, &info);

  free(work);
  free(iwork);
  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(7, rblapack_q, rblapack_m, rblapack_w, rblapack_z, rblapack_ifail, rblapack_info, rblapack_ab);
}

void
init_lapack_ssbevx(VALUE mLapack){
  rb_define_module_function(mLapack, "ssbevx", rblapack_ssbevx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
