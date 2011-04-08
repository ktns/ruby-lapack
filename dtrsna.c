#include "rb_lapack.h"

extern VOID dtrsna_(char *job, char *howmny, logical *select, integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, doublereal *s, doublereal *sep, integer *mm, integer *m, doublereal *work, integer *ldwork, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dtrsna(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  char job; 
  VALUE rblapack_howmny;
  char howmny; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_t;
  doublereal *t; 
  VALUE rblapack_vl;
  doublereal *vl; 
  VALUE rblapack_vr;
  doublereal *vr; 
  VALUE rblapack_ldwork;
  integer ldwork; 
  VALUE rblapack_s;
  doublereal *s; 
  VALUE rblapack_sep;
  doublereal *sep; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_info;
  integer info; 
  doublereal *work;
  integer *iwork;

  integer n;
  integer ldt;
  integer ldvl;
  integer ldvr;
  integer mm;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, sep, m, info = NumRu::Lapack.dtrsna( job, howmny, select, t, vl, vr, ldwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DTRSNA estimates reciprocal condition numbers for specified\n*  eigenvalues and/or right eigenvectors of a real upper\n*  quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q\n*  orthogonal).\n*\n*  T must be in Schur canonical form (as returned by DHSEQR), that is,\n*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each\n*  2-by-2 diagonal block has its diagonal elements equal and its\n*  off-diagonal elements of opposite sign.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies whether condition numbers are required for\n*          eigenvalues (S) or eigenvectors (SEP):\n*          = 'E': for eigenvalues only (S);\n*          = 'V': for eigenvectors only (SEP);\n*          = 'B': for both eigenvalues and eigenvectors (S and SEP).\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A': compute condition numbers for all eigenpairs;\n*          = 'S': compute condition numbers for selected eigenpairs\n*                 specified by the array SELECT.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which\n*          condition numbers are required. To select condition numbers\n*          for the eigenpair corresponding to a real eigenvalue w(j),\n*          SELECT(j) must be set to .TRUE.. To select condition numbers\n*          corresponding to a complex conjugate pair of eigenvalues w(j)\n*          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be\n*          set to .TRUE..\n*          If HOWMNY = 'A', SELECT is not referenced.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input) DOUBLE PRECISION array, dimension (LDT,N)\n*          The upper quasi-triangular matrix T, in Schur canonical form.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  VL      (input) DOUBLE PRECISION array, dimension (LDVL,M)\n*          If JOB = 'E' or 'B', VL must contain left eigenvectors of T\n*          (or of any Q*T*Q**T with Q orthogonal), corresponding to the\n*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors\n*          must be stored in consecutive columns of VL, as returned by\n*          DHSEIN or DTREVC.\n*          If JOB = 'V', VL is not referenced.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.\n*          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.\n*\n*  VR      (input) DOUBLE PRECISION array, dimension (LDVR,M)\n*          If JOB = 'E' or 'B', VR must contain right eigenvectors of T\n*          (or of any Q*T*Q**T with Q orthogonal), corresponding to the\n*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors\n*          must be stored in consecutive columns of VR, as returned by\n*          DHSEIN or DTREVC.\n*          If JOB = 'V', VR is not referenced.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.\n*          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.\n*\n*  S       (output) DOUBLE PRECISION array, dimension (MM)\n*          If JOB = 'E' or 'B', the reciprocal condition numbers of the\n*          selected eigenvalues, stored in consecutive elements of the\n*          array. For a complex conjugate pair of eigenvalues two\n*          consecutive elements of S are set to the same value. Thus\n*          S(j), SEP(j), and the j-th columns of VL and VR all\n*          correspond to the same eigenpair (but not in general the\n*          j-th eigenpair, unless all eigenpairs are selected).\n*          If JOB = 'V', S is not referenced.\n*\n*  SEP     (output) DOUBLE PRECISION array, dimension (MM)\n*          If JOB = 'V' or 'B', the estimated reciprocal condition\n*          numbers of the selected eigenvectors, stored in consecutive\n*          elements of the array. For a complex eigenvector two\n*          consecutive elements of SEP are set to the same value. If\n*          the eigenvalues cannot be reordered to compute SEP(j), SEP(j)\n*          is set to 0; this can only occur when the true value would be\n*          very small anyway.\n*          If JOB = 'E', SEP is not referenced.\n*\n*  MM      (input) INTEGER\n*          The number of elements in the arrays S (if JOB = 'E' or 'B')\n*           and/or SEP (if JOB = 'V' or 'B'). MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of elements of the arrays S and/or SEP actually\n*          used to store the estimated condition numbers.\n*          If HOWMNY = 'A', M is set to N.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,N+6)\n*          If JOB = 'E', WORK is not referenced.\n*\n*  LDWORK  (input) INTEGER\n*          The leading dimension of the array WORK.\n*          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.\n*\n*  IWORK   (workspace) INTEGER array, dimension (2*(N-1))\n*          If JOB = 'E', IWORK is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The reciprocal of the condition number of an eigenvalue lambda is\n*  defined as\n*\n*          S(lambda) = |v'*u| / (norm(u)*norm(v))\n*\n*  where u and v are the right and left eigenvectors of T corresponding\n*  to lambda; v' denotes the conjugate-transpose of v, and norm(u)\n*  denotes the Euclidean norm. These reciprocal condition numbers always\n*  lie between zero (very badly conditioned) and one (very well\n*  conditioned). If n = 1, S(lambda) is defined to be 1.\n*\n*  An approximate error bound for a computed eigenvalue W(i) is given by\n*\n*                      EPS * norm(T) / S(i)\n*\n*  where EPS is the machine precision.\n*\n*  The reciprocal of the condition number of the right eigenvector u\n*  corresponding to lambda is defined as follows. Suppose\n*\n*              T = ( lambda  c  )\n*                  (   0    T22 )\n*\n*  Then the reciprocal condition number is\n*\n*          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )\n*\n*  where sigma-min denotes the smallest singular value. We approximate\n*  the smallest singular value by the reciprocal of an estimate of the\n*  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is\n*  defined to be abs(T(1,1)).\n*\n*  An approximate error bound for a computed right eigenvector VR(i)\n*  is given by\n*\n*                      EPS * norm(T) / SEP(i)\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, sep, m, info = NumRu::Lapack.dtrsna( job, howmny, select, t, vl, vr, ldwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_job = argv[0];
  rblapack_howmny = argv[1];
  rblapack_select = argv[2];
  rblapack_t = argv[3];
  rblapack_vl = argv[4];
  rblapack_vr = argv[5];
  rblapack_ldwork = argv[6];
  if (rb_options != Qnil) {
  }

  howmny = StringValueCStr(rblapack_howmny)[0];
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (4th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_t);
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_DFLOAT)
    rblapack_t = na_change_type(rblapack_t, NA_DFLOAT);
  t = NA_PTR_TYPE(rblapack_t, doublereal*);
  if (!NA_IsNArray(rblapack_vl))
    rb_raise(rb_eArgError, "vl (5th argument) must be NArray");
  if (NA_RANK(rblapack_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (5th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_vl);
  ldvl = NA_SHAPE0(rblapack_vl);
  if (NA_TYPE(rblapack_vl) != NA_DFLOAT)
    rblapack_vl = na_change_type(rblapack_vl, NA_DFLOAT);
  vl = NA_PTR_TYPE(rblapack_vl, doublereal*);
  if (!NA_IsNArray(rblapack_vr))
    rb_raise(rb_eArgError, "vr (6th argument) must be NArray");
  if (NA_RANK(rblapack_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_vr) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rblapack_vr);
  if (NA_TYPE(rblapack_vr) != NA_DFLOAT)
    rblapack_vr = na_change_type(rblapack_vr, NA_DFLOAT);
  vr = NA_PTR_TYPE(rblapack_vr, doublereal*);
  if (!NA_IsNArray(rblapack_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rblapack_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of t");
  if (NA_TYPE(rblapack_select) != NA_LINT)
    rblapack_select = na_change_type(rblapack_select, NA_LINT);
  select = NA_PTR_TYPE(rblapack_select, logical*);
  job = StringValueCStr(rblapack_job)[0];
  ldwork = ((lsame_(&job,"V")) || (lsame_(&job,"B"))) ? n : 1;
  mm = m;
  {
    int shape[1];
    shape[0] = mm;
    rblapack_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rblapack_s, doublereal*);
  {
    int shape[1];
    shape[0] = mm;
    rblapack_sep = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sep = NA_PTR_TYPE(rblapack_sep, doublereal*);
  work = ALLOC_N(doublereal, (lsame_(&job,"E") ? 0 : ldwork)*(lsame_(&job,"E") ? 0 : n+6));
  iwork = ALLOC_N(integer, (lsame_(&job,"E") ? 0 : 2*(n-1)));

  dtrsna_(&job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, &m, work, &ldwork, iwork, &info);

  free(work);
  free(iwork);
  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_s, rblapack_sep, rblapack_m, rblapack_info);
}

void
init_lapack_dtrsna(VALUE mLapack){
  rb_define_module_function(mLapack, "dtrsna", rblapack_dtrsna, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
