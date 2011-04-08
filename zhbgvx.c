#include "rb_lapack.h"

extern VOID zhbgvx_(char *jobz, char *range, char *uplo, integer *n, integer *ka, integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, integer *ldbb, doublecomplex *q, integer *ldq, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z, integer *ldz, doublecomplex *work, doublereal *rwork, integer *iwork, integer *ifail, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zhbgvx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobz;
  char jobz; 
  VALUE rblapack_range;
  char range; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ka;
  integer ka; 
  VALUE rblapack_kb;
  integer kb; 
  VALUE rblapack_ab;
  doublecomplex *ab; 
  VALUE rblapack_bb;
  doublecomplex *bb; 
  VALUE rblapack_vl;
  doublereal vl; 
  VALUE rblapack_vu;
  doublereal vu; 
  VALUE rblapack_il;
  integer il; 
  VALUE rblapack_iu;
  integer iu; 
  VALUE rblapack_abstol;
  doublereal abstol; 
  VALUE rblapack_q;
  doublecomplex *q; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_w;
  doublereal *w; 
  VALUE rblapack_z;
  doublecomplex *z; 
  VALUE rblapack_ifail;
  integer *ifail; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ab_out__;
  doublecomplex *ab_out__;
  VALUE rblapack_bb_out__;
  doublecomplex *bb_out__;
  doublecomplex *work;
  doublereal *rwork;
  integer *iwork;

  integer ldab;
  integer n;
  integer ldbb;
  integer ldq;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  q, m, w, z, ifail, info, ab, bb = NumRu::Lapack.zhbgvx( jobz, range, uplo, ka, kb, ab, bb, vl, vu, il, iu, abstol, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZHBGVX computes all the eigenvalues, and optionally, the eigenvectors\n*  of a complex generalized Hermitian-definite banded eigenproblem, of\n*  the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian\n*  and banded, and B is also positive definite.  Eigenvalues and\n*  eigenvectors can be selected by specifying either all eigenvalues,\n*  a range of values or a range of indices for the desired eigenvalues.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found;\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found;\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangles of A and B are stored;\n*          = 'L':  Lower triangles of A and B are stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  KA      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'. KA >= 0.\n*\n*  KB      (input) INTEGER\n*          The number of superdiagonals of the matrix B if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'. KB >= 0.\n*\n*  AB      (input/output) COMPLEX*16 array, dimension (LDAB, N)\n*          On entry, the upper or lower triangle of the Hermitian band\n*          matrix A, stored in the first ka+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).\n*\n*          On exit, the contents of AB are destroyed.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KA+1.\n*\n*  BB      (input/output) COMPLEX*16 array, dimension (LDBB, N)\n*          On entry, the upper or lower triangle of the Hermitian band\n*          matrix B, stored in the first kb+1 rows of the array.  The\n*          j-th column of B is stored in the j-th column of the array BB\n*          as follows:\n*          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;\n*          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).\n*\n*          On exit, the factor S from the split Cholesky factorization\n*          B = S**H*S, as returned by ZPBSTF.\n*\n*  LDBB    (input) INTEGER\n*          The leading dimension of the array BB.  LDBB >= KB+1.\n*\n*  Q       (output) COMPLEX*16 array, dimension (LDQ, N)\n*          If JOBZ = 'V', the n-by-n matrix used in the reduction of\n*          A*x = (lambda)*B*x to standard form, i.e. C*x = (lambda)*x,\n*          and consequently C to tridiagonal form.\n*          If JOBZ = 'N', the array Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  If JOBZ = 'N',\n*          LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) DOUBLE PRECISION\n*          The absolute error tolerance for the eigenvalues.\n*          An approximate eigenvalue is accepted as converged\n*          when it is determined to lie in an interval [a,b]\n*          of width less than or equal to\n*\n*                  ABSTOL + EPS *   max( |a|,|b| ) ,\n*\n*          where EPS is the machine precision.  If ABSTOL is less than\n*          or equal to zero, then  EPS*|T|  will be used in its place,\n*          where |T| is the 1-norm of the tridiagonal matrix obtained\n*          by reducing AP to tridiagonal form.\n*\n*          Eigenvalues will be computed most accurately when ABSTOL is\n*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.\n*          If this routine returns with INFO>0, indicating that some\n*          eigenvectors did not converge, try setting ABSTOL to\n*          2*DLAMCH('S').\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of\n*          eigenvectors, with the i-th column of Z holding the\n*          eigenvector associated with W(i). The eigenvectors are\n*          normalized so that Z**H*B*Z = I.\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= N.\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (7*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (5*N)\n*\n*  IFAIL   (output) INTEGER array, dimension (N)\n*          If JOBZ = 'V', then if INFO = 0, the first M elements of\n*          IFAIL are zero.  If INFO > 0, then IFAIL contains the\n*          indices of the eigenvectors that failed to converge.\n*          If JOBZ = 'N', then IFAIL is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is:\n*             <= N:  then i eigenvectors failed to converge.  Their\n*                    indices are stored in array IFAIL.\n*             > N:   if INFO = N + i, for 1 <= i <= N, then ZPBSTF\n*                    returned INFO = i: B is not positive definite.\n*                    The factorization of B could not be completed and\n*                    no eigenvalues or eigenvectors were computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  q, m, w, z, ifail, info, ab, bb = NumRu::Lapack.zhbgvx( jobz, range, uplo, ka, kb, ab, bb, vl, vu, il, iu, abstol, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rblapack_jobz = argv[0];
  rblapack_range = argv[1];
  rblapack_uplo = argv[2];
  rblapack_ka = argv[3];
  rblapack_kb = argv[4];
  rblapack_ab = argv[5];
  rblapack_bb = argv[6];
  rblapack_vl = argv[7];
  rblapack_vu = argv[8];
  rblapack_il = argv[9];
  rblapack_iu = argv[10];
  rblapack_abstol = argv[11];
  if (rb_options != Qnil) {
  }

  abstol = NUM2DBL(rblapack_abstol);
  vl = NUM2DBL(rblapack_vl);
  if (!NA_IsNArray(rblapack_bb))
    rb_raise(rb_eArgError, "bb (7th argument) must be NArray");
  if (NA_RANK(rblapack_bb) != 2)
    rb_raise(rb_eArgError, "rank of bb (7th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_bb);
  ldbb = NA_SHAPE0(rblapack_bb);
  if (NA_TYPE(rblapack_bb) != NA_DCOMPLEX)
    rblapack_bb = na_change_type(rblapack_bb, NA_DCOMPLEX);
  bb = NA_PTR_TYPE(rblapack_bb, doublecomplex*);
  ka = NUM2INT(rblapack_ka);
  iu = NUM2INT(rblapack_iu);
  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (6th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 1 of bb");
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_DCOMPLEX)
    rblapack_ab = na_change_type(rblapack_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rblapack_ab, doublecomplex*);
  kb = NUM2INT(rblapack_kb);
  vu = NUM2DBL(rblapack_vu);
  jobz = StringValueCStr(rblapack_jobz)[0];
  il = NUM2INT(rblapack_il);
  range = StringValueCStr(rblapack_range)[0];
  uplo = StringValueCStr(rblapack_uplo)[0];
  ldq = 1 ? jobz = 'n' : max(1,n) ? jobz = 'v' : 0;
  ldz = lsame_(&jobz,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rblapack_q, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, doublereal*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, doublecomplex*);
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
    rblapack_ab_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rblapack_ab_out__, doublecomplex*);
  MEMCPY(ab_out__, ab, doublecomplex, NA_TOTAL(rblapack_ab));
  rblapack_ab = rblapack_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldbb;
    shape[1] = n;
    rblapack_bb_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  bb_out__ = NA_PTR_TYPE(rblapack_bb_out__, doublecomplex*);
  MEMCPY(bb_out__, bb, doublecomplex, NA_TOTAL(rblapack_bb));
  rblapack_bb = rblapack_bb_out__;
  bb = bb_out__;
  work = ALLOC_N(doublecomplex, (n));
  rwork = ALLOC_N(doublereal, (7*n));
  iwork = ALLOC_N(integer, (5*n));

  zhbgvx_(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, rwork, iwork, ifail, &info);

  free(work);
  free(rwork);
  free(iwork);
  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(8, rblapack_q, rblapack_m, rblapack_w, rblapack_z, rblapack_ifail, rblapack_info, rblapack_ab, rblapack_bb);
}

void
init_lapack_zhbgvx(VALUE mLapack){
  rb_define_module_function(mLapack, "zhbgvx", rblapack_zhbgvx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
