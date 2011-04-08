#include "rb_lapack.h"

extern VOID dspevd_(char *jobz, char *uplo, integer *n, doublereal *ap, doublereal *w, doublereal *z, integer *ldz, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dspevd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobz;
  char jobz; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ap;
  doublereal *ap; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_liwork;
  integer liwork; 
  VALUE rblapack_w;
  doublereal *w; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ap_out__;
  doublereal *ap_out__;

  integer ldap;
  integer n;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, z, work, iwork, info, ap = NumRu::Lapack.dspevd( jobz, uplo, ap, lwork, liwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSPEVD computes all the eigenvalues and, optionally, eigenvectors\n*  of a real symmetric matrix A in packed storage. If eigenvectors are\n*  desired, it uses a divide and conquer algorithm.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the symmetric matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*\n*          On exit, AP is overwritten by values generated during the\n*          reduction to tridiagonal form.  If UPLO = 'U', the diagonal\n*          and first superdiagonal of the tridiagonal matrix T overwrite\n*          the corresponding elements of A, and if UPLO = 'L', the\n*          diagonal and first subdiagonal of T overwrite the\n*          corresponding elements of A.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal\n*          eigenvectors of the matrix A, with the i-th column of Z\n*          holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array,\n*                                         dimension (LWORK)\n*          On exit, if INFO = 0, WORK(1) returns the required LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If N <= 1,               LWORK must be at least 1.\n*          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N.\n*          If JOBZ = 'V' and N > 1, LWORK must be at least\n*                                                 1 + 6*N + N**2.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the required sizes of the WORK and IWORK\n*          arrays, returns these values as the first entries of the WORK\n*          and IWORK arrays, and no error message related to LWORK or\n*          LIWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the required LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.\n*          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.\n*          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the required sizes of the WORK and\n*          IWORK arrays, returns these values as the first entries of\n*          the WORK and IWORK arrays, and no error message related to\n*          LWORK or LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, the algorithm failed to converge; i\n*                off-diagonal elements of an intermediate tridiagonal\n*                form did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, z, work, iwork, info, ap = NumRu::Lapack.dspevd( jobz, uplo, ap, lwork, liwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_jobz = argv[0];
  rblapack_uplo = argv[1];
  rblapack_ap = argv[2];
  rblapack_lwork = argv[3];
  rblapack_liwork = argv[4];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  jobz = StringValueCStr(rblapack_jobz)[0];
  liwork = NUM2INT(rblapack_liwork);
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  ldap = NA_SHAPE0(rblapack_ap);
  if (NA_TYPE(rblapack_ap) != NA_DFLOAT)
    rblapack_ap = na_change_type(rblapack_ap, NA_DFLOAT);
  ap = NA_PTR_TYPE(rblapack_ap, doublereal*);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
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
    rblapack_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rblapack_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
  {
    int shape[1];
    shape[0] = ldap;
    rblapack_ap_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rblapack_ap_out__, doublereal*);
  MEMCPY(ap_out__, ap, doublereal, NA_TOTAL(rblapack_ap));
  rblapack_ap = rblapack_ap_out__;
  ap = ap_out__;

  dspevd_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, iwork, &liwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_w, rblapack_z, rblapack_work, rblapack_iwork, rblapack_info, rblapack_ap);
}

void
init_lapack_dspevd(VALUE mLapack){
  rb_define_module_function(mLapack, "dspevd", rblapack_dspevd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
