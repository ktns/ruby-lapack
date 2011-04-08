#include "rb_lapack.h"

static logical
rblapack_select(real *arg0, real *arg1){
  VALUE rblapack_arg0, rblapack_arg1;

  VALUE rblapack_ret;
  logical ret;

  rblapack_arg0 = rb_float_new((double)(*arg0));
  rblapack_arg1 = rb_float_new((double)(*arg1));

  rblapack_ret = rb_yield_values(2, rblapack_arg0, rblapack_arg1);

  ret = (rblapack_ret == Qtrue);
  return ret;
}

extern VOID sgeesx_(char *jobvs, char *sort, L_fp *select, char *sense, integer *n, real *a, integer *lda, integer *sdim, real *wr, real *wi, real *vs, integer *ldvs, real *rconde, real *rcondv, real *work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sgeesx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobvs;
  char jobvs; 
  VALUE rblapack_sort;
  char sort; 
  VALUE rblapack_sense;
  char sense; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_liwork;
  integer liwork; 
  VALUE rblapack_sdim;
  integer sdim; 
  VALUE rblapack_wr;
  real *wr; 
  VALUE rblapack_wi;
  real *wi; 
  VALUE rblapack_vs;
  real *vs; 
  VALUE rblapack_rconde;
  real rconde; 
  VALUE rblapack_rcondv;
  real rcondv; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  real *a_out__;
  logical *bwork;

  integer lda;
  integer n;
  integer ldvs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sdim, wr, wi, vs, rconde, rcondv, work, iwork, info, a = NumRu::Lapack.sgeesx( jobvs, sort, sense, a, lwork, liwork, [:usage => usage, :help => help]){|a,b| ... }\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGEESX computes for an N-by-N real nonsymmetric matrix A, the\n*  eigenvalues, the real Schur form T, and, optionally, the matrix of\n*  Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).\n*\n*  Optionally, it also orders the eigenvalues on the diagonal of the\n*  real Schur form so that selected eigenvalues are at the top left;\n*  computes a reciprocal condition number for the average of the\n*  selected eigenvalues (RCONDE); and computes a reciprocal condition\n*  number for the right invariant subspace corresponding to the\n*  selected eigenvalues (RCONDV).  The leading columns of Z form an\n*  orthonormal basis for this invariant subspace.\n*\n*  For further explanation of the reciprocal condition numbers RCONDE\n*  and RCONDV, see Section 4.10 of the LAPACK Users' Guide (where\n*  these quantities are called s and sep respectively).\n*\n*  A real matrix is in real Schur form if it is upper quasi-triangular\n*  with 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in\n*  the form\n*            [  a  b  ]\n*            [  c  a  ]\n*\n*  where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVS   (input) CHARACTER*1\n*          = 'N': Schur vectors are not computed;\n*          = 'V': Schur vectors are computed.\n*\n*  SORT    (input) CHARACTER*1\n*          Specifies whether or not to order the eigenvalues on the\n*          diagonal of the Schur form.\n*          = 'N': Eigenvalues are not ordered;\n*          = 'S': Eigenvalues are ordered (see SELECT).\n*\n*  SELECT  (external procedure) LOGICAL FUNCTION of two REAL arguments\n*          SELECT must be declared EXTERNAL in the calling subroutine.\n*          If SORT = 'S', SELECT is used to select eigenvalues to sort\n*          to the top left of the Schur form.\n*          If SORT = 'N', SELECT is not referenced.\n*          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if\n*          SELECT(WR(j),WI(j)) is true; i.e., if either one of a\n*          complex conjugate pair of eigenvalues is selected, then both\n*          are.  Note that a selected complex eigenvalue may no longer\n*          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since\n*          ordering may change the value of complex eigenvalues\n*          (especially if the eigenvalue is ill-conditioned); in this\n*          case INFO may be set to N+3 (see INFO below).\n*\n*  SENSE   (input) CHARACTER*1\n*          Determines which reciprocal condition numbers are computed.\n*          = 'N': None are computed;\n*          = 'E': Computed for average of selected eigenvalues only;\n*          = 'V': Computed for selected right invariant subspace only;\n*          = 'B': Computed for both.\n*          If SENSE = 'E', 'V' or 'B', SORT must equal 'S'.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA, N)\n*          On entry, the N-by-N matrix A.\n*          On exit, A is overwritten by its real Schur form T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  SDIM    (output) INTEGER\n*          If SORT = 'N', SDIM = 0.\n*          If SORT = 'S', SDIM = number of eigenvalues (after sorting)\n*                         for which SELECT is true. (Complex conjugate\n*                         pairs for which SELECT is true for either\n*                         eigenvalue count as 2.)\n*\n*  WR      (output) REAL array, dimension (N)\n*  WI      (output) REAL array, dimension (N)\n*          WR and WI contain the real and imaginary parts, respectively,\n*          of the computed eigenvalues, in the same order that they\n*          appear on the diagonal of the output Schur form T.  Complex\n*          conjugate pairs of eigenvalues appear consecutively with the\n*          eigenvalue having the positive imaginary part first.\n*\n*  VS      (output) REAL array, dimension (LDVS,N)\n*          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur\n*          vectors.\n*          If JOBVS = 'N', VS is not referenced.\n*\n*  LDVS    (input) INTEGER\n*          The leading dimension of the array VS.  LDVS >= 1, and if\n*          JOBVS = 'V', LDVS >= N.\n*\n*  RCONDE  (output) REAL\n*          If SENSE = 'E' or 'B', RCONDE contains the reciprocal\n*          condition number for the average of the selected eigenvalues.\n*          Not referenced if SENSE = 'N' or 'V'.\n*\n*  RCONDV  (output) REAL\n*          If SENSE = 'V' or 'B', RCONDV contains the reciprocal\n*          condition number for the selected right invariant subspace.\n*          Not referenced if SENSE = 'N' or 'E'.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,3*N).\n*          Also, if SENSE = 'E' or 'V' or 'B',\n*          LWORK >= N+2*SDIM*(N-SDIM), where SDIM is the number of\n*          selected eigenvalues computed by this routine.  Note that\n*          N+2*SDIM*(N-SDIM) <= N+N*N/2. Note also that an error is only\n*          returned if LWORK < max(1,3*N), but if SENSE = 'E' or 'V' or\n*          'B' this may not be large enough.\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates upper bounds on the optimal sizes of the\n*          arrays WORK and IWORK, returns these values as the first\n*          entries of the WORK and IWORK arrays, and no error messages\n*          related to LWORK or LIWORK are issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.\n*          LIWORK >= 1; if SENSE = 'V' or 'B', LIWORK >= SDIM*(N-SDIM).\n*          Note that SDIM*(N-SDIM) <= N*N/4. Note also that an error is\n*          only returned if LIWORK < 1, but if SENSE = 'V' or 'B' this\n*          may not be large enough.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates upper bounds on the optimal sizes of\n*          the arrays WORK and IWORK, returns these values as the first\n*          entries of the WORK and IWORK arrays, and no error messages\n*          related to LWORK or LIWORK are issued by XERBLA.\n*\n*  BWORK   (workspace) LOGICAL array, dimension (N)\n*          Not referenced if SORT = 'N'.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value.\n*          > 0: if INFO = i, and i is\n*             <= N: the QR algorithm failed to compute all the\n*                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI\n*                   contain those eigenvalues which have converged; if\n*                   JOBVS = 'V', VS contains the transformation which\n*                   reduces A to its partially converged Schur form.\n*             = N+1: the eigenvalues could not be reordered because some\n*                   eigenvalues were too close to separate (the problem\n*                   is very ill-conditioned);\n*             = N+2: after reordering, roundoff changed values of some\n*                   complex eigenvalues so that leading eigenvalues in\n*                   the Schur form no longer satisfy SELECT=.TRUE.  This\n*                   could also be caused by underflow due to scaling.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sdim, wr, wi, vs, rconde, rcondv, work, iwork, info, a = NumRu::Lapack.sgeesx( jobvs, sort, sense, a, lwork, liwork, [:usage => usage, :help => help]){|a,b| ... }\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_jobvs = argv[0];
  rblapack_sort = argv[1];
  rblapack_sense = argv[2];
  rblapack_a = argv[3];
  rblapack_lwork = argv[4];
  rblapack_liwork = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  jobvs = StringValueCStr(rblapack_jobvs)[0];
  sort = StringValueCStr(rblapack_sort)[0];
  lwork = NUM2INT(rblapack_lwork);
  sense = StringValueCStr(rblapack_sense)[0];
  liwork = NUM2INT(rblapack_liwork);
  ldvs = lsame_(&jobvs,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_wr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rblapack_wr, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_wi = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rblapack_wi, real*);
  {
    int shape[2];
    shape[0] = ldvs;
    shape[1] = n;
    rblapack_vs = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vs = NA_PTR_TYPE(rblapack_vs, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rblapack_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  sgeesx_(&jobvs, &sort, rblapack_select, &sense, &n, a, &lda, &sdim, wr, wi, vs, &ldvs, &rconde, &rcondv, work, &lwork, iwork, &liwork, bwork, &info);

  free(bwork);
  rblapack_sdim = INT2NUM(sdim);
  rblapack_rconde = rb_float_new((double)rconde);
  rblapack_rcondv = rb_float_new((double)rcondv);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(10, rblapack_sdim, rblapack_wr, rblapack_wi, rblapack_vs, rblapack_rconde, rblapack_rcondv, rblapack_work, rblapack_iwork, rblapack_info, rblapack_a);
}

void
init_lapack_sgeesx(VALUE mLapack){
  rb_define_module_function(mLapack, "sgeesx", rblapack_sgeesx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
