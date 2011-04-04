#include "rb_lapack.h"

extern VOID sgsvj0_(char *jobv, integer *m, integer *n, real *a, integer *lda, real *d, real *sva, integer *mv, real *v, integer *ldv, integer *eps, integer *sfmin, real *tol, integer *nsweep, real *work, integer *lwork, integer *info);

static VALUE
rb_sgsvj0(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  real *a; 
  VALUE rb_d;
  real *d; 
  VALUE rb_sva;
  real *sva; 
  VALUE rb_mv;
  integer mv; 
  VALUE rb_v;
  real *v; 
  VALUE rb_eps;
  integer eps; 
  VALUE rb_sfmin;
  integer sfmin; 
  VALUE rb_tol;
  real tol; 
  VALUE rb_nsweep;
  integer nsweep; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_sva_out__;
  real *sva_out__;
  VALUE rb_v_out__;
  real *v_out__;
  real *work;

  integer lda;
  integer n;
  integer ldv;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a, d, sva, v = NumRu::Lapack.sgsvj0( jobv, m, a, d, sva, mv, v, eps, sfmin, tol, nsweep)\n    or\n  NumRu::Lapack.sgsvj0  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGSVJ0 is called from SGESVJ as a pre-processor and that is its main\n*  purpose. It applies Jacobi rotations in the same way as SGESVJ does, but\n*  it does not check convergence (stopping criterion). Few tuning\n*  parameters (marked by [TP]) are available for the implementer.\n*\n*  Further Details\n*  ~~~~~~~~~~~~~~~\n*  SGSVJ0 is used just to enable SGESVJ to call a simplified version of\n*  itself to work on a submatrix of the original matrix.\n*\n*  Contributors\n*  ~~~~~~~~~~~~\n*  Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)\n*\n*  Bugs, Examples and Comments\n*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~\n*  Please report all bugs and send interesting test examples and comments to\n*  drmac@math.hr. Thank you.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBV    (input) CHARACTER*1\n*          Specifies whether the output from this procedure is used\n*          to compute the matrix V:\n*          = 'V': the product of the Jacobi rotations is accumulated\n*                 by postmulyiplying the N-by-N array V.\n*                (See the description of V.)\n*          = 'A': the product of the Jacobi rotations is accumulated\n*                 by postmulyiplying the MV-by-N array V.\n*                (See the descriptions of MV and V.)\n*          = 'N': the Jacobi rotations are not accumulated.\n*\n*  M       (input) INTEGER\n*          The number of rows of the input matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the input matrix A.\n*          M >= N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, M-by-N matrix A, such that A*diag(D) represents\n*          the input matrix.\n*          On exit,\n*          A_onexit * D_onexit represents the input matrix A*diag(D)\n*          post-multiplied by a sequence of Jacobi rotations, where the\n*          rotation threshold and the total number of sweeps are given in\n*          TOL and NSWEEP, respectively.\n*          (See the descriptions of D, TOL and NSWEEP.)\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  D       (input/workspace/output) REAL array, dimension (N)\n*          The array D accumulates the scaling factors from the fast scaled\n*          Jacobi rotations.\n*          On entry, A*diag(D) represents the input matrix.\n*          On exit, A_onexit*diag(D_onexit) represents the input matrix\n*          post-multiplied by a sequence of Jacobi rotations, where the\n*          rotation threshold and the total number of sweeps are given in\n*          TOL and NSWEEP, respectively.\n*          (See the descriptions of A, TOL and NSWEEP.)\n*\n*  SVA     (input/workspace/output) REAL array, dimension (N)\n*          On entry, SVA contains the Euclidean norms of the columns of\n*          the matrix A*diag(D).\n*          On exit, SVA contains the Euclidean norms of the columns of\n*          the matrix onexit*diag(D_onexit).\n*\n*  MV      (input) INTEGER\n*          If JOBV .EQ. 'A', then MV rows of V are post-multipled by a\n*                           sequence of Jacobi rotations.\n*          If JOBV = 'N',   then MV is not referenced.\n*\n*  V       (input/output) REAL array, dimension (LDV,N)\n*          If JOBV .EQ. 'V' then N rows of V are post-multipled by a\n*                           sequence of Jacobi rotations.\n*          If JOBV .EQ. 'A' then MV rows of V are post-multipled by a\n*                           sequence of Jacobi rotations.\n*          If JOBV = 'N',   then V is not referenced.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V,  LDV >= 1.\n*          If JOBV = 'V', LDV .GE. N.\n*          If JOBV = 'A', LDV .GE. MV.\n*\n*  EPS     (input) INTEGER\n*          EPS = SLAMCH('Epsilon')\n*\n*  SFMIN   (input) INTEGER\n*          SFMIN = SLAMCH('Safe Minimum')\n*\n*  TOL     (input) REAL\n*          TOL is the threshold for Jacobi rotations. For a pair\n*          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is\n*          applied only if ABS(COS(angle(A(:,p),A(:,q)))) .GT. TOL.\n*\n*  NSWEEP  (input) INTEGER\n*          NSWEEP is the number of sweeps of Jacobi rotations to be\n*          performed.\n*\n*  WORK    (workspace) REAL array, dimension LWORK.\n*\n*  LWORK   (input) INTEGER\n*          LWORK is the dimension of WORK. LWORK .GE. M.\n*\n*  INFO    (output) INTEGER\n*          = 0 : successful exit.\n*          < 0 : if INFO = -i, then the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n*     .. Local Parameters ..\n      REAL               ZERO, HALF, ONE, TWO\n      PARAMETER          ( ZERO = 0.0E0, HALF = 0.5E0, ONE = 1.0E0,\n     +                   TWO = 2.0E0 )\n*     ..\n*     .. Local Scalars ..\n      REAL               AAPP, AAPP0, AAPQ, AAQQ, APOAQ, AQOAP, BIG,\n     +                   BIGTHETA, CS, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS,\n     +                   ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA,\n     +                   THSIGN\n      INTEGER            BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1,\n     +                   ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, NBL,\n     +                   NOTROT, p, PSKIPPED, q, ROWSKIP, SWBAND\n      LOGICAL            APPLV, ROTOK, RSVEC\n*     ..\n*     .. Local Arrays ..\n      REAL               FASTR( 5 )\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, AMAX1, AMIN1, FLOAT, MIN0, SIGN, SQRT\n*     ..\n*     .. External Functions ..\n      REAL               SDOT, SNRM2\n      INTEGER            ISAMAX\n      LOGICAL            LSAME\n      EXTERNAL           ISAMAX, LSAME, SDOT, SNRM2\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SAXPY, SCOPY, SLASCL, SLASSQ, SROTM, SSWAP\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_jobv = argv[0];
  rb_m = argv[1];
  rb_a = argv[2];
  rb_d = argv[3];
  rb_sva = argv[4];
  rb_mv = argv[5];
  rb_v = argv[6];
  rb_eps = argv[7];
  rb_sfmin = argv[8];
  rb_tol = argv[9];
  rb_nsweep = argv[10];

  sfmin = NUM2INT(rb_sfmin);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (7th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_v);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_SFLOAT)
    rb_v = na_change_type(rb_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rb_v, real*);
  nsweep = NUM2INT(rb_nsweep);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 1 of v");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of v");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_sva))
    rb_raise(rb_eArgError, "sva (5th argument) must be NArray");
  if (NA_RANK(rb_sva) != 1)
    rb_raise(rb_eArgError, "rank of sva (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_sva) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of sva must be the same as shape 1 of v");
  if (NA_TYPE(rb_sva) != NA_SFLOAT)
    rb_sva = na_change_type(rb_sva, NA_SFLOAT);
  sva = NA_PTR_TYPE(rb_sva, real*);
  jobv = StringValueCStr(rb_jobv)[0];
  tol = (real)NUM2DBL(rb_tol);
  eps = NUM2INT(rb_eps);
  mv = NUM2INT(rb_mv);
  lwork = m;
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_sva_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  sva_out__ = NA_PTR_TYPE(rb_sva_out__, real*);
  MEMCPY(sva_out__, sva, real, NA_TOTAL(rb_sva));
  rb_sva = rb_sva_out__;
  sva = sva_out__;
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = n;
    rb_v_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, real*);
  MEMCPY(v_out__, v, real, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;
  work = ALLOC_N(real, (lwork));

  sgsvj0_(&jobv, &m, &n, a, &lda, d, sva, &mv, v, &ldv, &eps, &sfmin, &tol, &nsweep, work, &lwork, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_info, rb_a, rb_d, rb_sva, rb_v);
}

void
init_lapack_sgsvj0(VALUE mLapack){
  rb_define_module_function(mLapack, "sgsvj0", rb_sgsvj0, -1);
}
