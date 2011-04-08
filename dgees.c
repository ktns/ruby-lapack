#include "rb_lapack.h"

static logical
rblapack_select(doublereal *arg0, doublereal *arg1){
  VALUE rblapack_arg0, rblapack_arg1;

  VALUE rblapack_ret;
  logical ret;

  rblapack_arg0 = rb_float_new((double)(*arg0));
  rblapack_arg1 = rb_float_new((double)(*arg1));

  rblapack_ret = rb_yield_values(2, rblapack_arg0, rblapack_arg1);

  ret = (rblapack_ret == Qtrue);
  return ret;
}

extern VOID dgees_(char *jobvs, char *sort, L_fp *select, integer *n, doublereal *a, integer *lda, integer *sdim, doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work, integer *lwork, logical *bwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgees(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobvs;
  char jobvs; 
  VALUE rblapack_sort;
  char sort; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_sdim;
  integer sdim; 
  VALUE rblapack_wr;
  doublereal *wr; 
  VALUE rblapack_wi;
  doublereal *wi; 
  VALUE rblapack_vs;
  doublereal *vs; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  logical *bwork;

  integer lda;
  integer n;
  integer ldvs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sdim, wr, wi, vs, work, info, a = NumRu::Lapack.dgees( jobvs, sort, a, lwork, [:usage => usage, :help => help]){|a,b| ... }\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI, VS, LDVS, WORK, LWORK, BWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGEES computes for an N-by-N real nonsymmetric matrix A, the\n*  eigenvalues, the real Schur form T, and, optionally, the matrix of\n*  Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).\n*\n*  Optionally, it also orders the eigenvalues on the diagonal of the\n*  real Schur form so that selected eigenvalues are at the top left.\n*  The leading columns of Z then form an orthonormal basis for the\n*  invariant subspace corresponding to the selected eigenvalues.\n*\n*  A matrix is in real Schur form if it is upper quasi-triangular with\n*  1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the\n*  form\n*          [  a  b  ]\n*          [  c  a  ]\n*\n*  where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVS   (input) CHARACTER*1\n*          = 'N': Schur vectors are not computed;\n*          = 'V': Schur vectors are computed.\n*\n*  SORT    (input) CHARACTER*1\n*          Specifies whether or not to order the eigenvalues on the\n*          diagonal of the Schur form.\n*          = 'N': Eigenvalues are not ordered;\n*          = 'S': Eigenvalues are ordered (see SELECT).\n*\n*  SELECT  (external procedure) LOGICAL FUNCTION of two DOUBLE PRECISION arguments\n*          SELECT must be declared EXTERNAL in the calling subroutine.\n*          If SORT = 'S', SELECT is used to select eigenvalues to sort\n*          to the top left of the Schur form.\n*          If SORT = 'N', SELECT is not referenced.\n*          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if\n*          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex\n*          conjugate pair of eigenvalues is selected, then both complex\n*          eigenvalues are selected.\n*          Note that a selected complex eigenvalue may no longer\n*          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since\n*          ordering may change the value of complex eigenvalues\n*          (especially if the eigenvalue is ill-conditioned); in this\n*          case INFO is set to N+2 (see INFO below).\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the N-by-N matrix A.\n*          On exit, A has been overwritten by its real Schur form T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  SDIM    (output) INTEGER\n*          If SORT = 'N', SDIM = 0.\n*          If SORT = 'S', SDIM = number of eigenvalues (after sorting)\n*                         for which SELECT is true. (Complex conjugate\n*                         pairs for which SELECT is true for either\n*                         eigenvalue count as 2.)\n*\n*  WR      (output) DOUBLE PRECISION array, dimension (N)\n*  WI      (output) DOUBLE PRECISION array, dimension (N)\n*          WR and WI contain the real and imaginary parts,\n*          respectively, of the computed eigenvalues in the same order\n*          that they appear on the diagonal of the output Schur form T.\n*          Complex conjugate pairs of eigenvalues will appear\n*          consecutively with the eigenvalue having the positive\n*          imaginary part first.\n*\n*  VS      (output) DOUBLE PRECISION array, dimension (LDVS,N)\n*          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur\n*          vectors.\n*          If JOBVS = 'N', VS is not referenced.\n*\n*  LDVS    (input) INTEGER\n*          The leading dimension of the array VS.  LDVS >= 1; if\n*          JOBVS = 'V', LDVS >= N.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) contains the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,3*N).\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  BWORK   (workspace) LOGICAL array, dimension (N)\n*          Not referenced if SORT = 'N'.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value.\n*          > 0: if INFO = i, and i is\n*             <= N: the QR algorithm failed to compute all the\n*                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI\n*                   contain those eigenvalues which have converged; if\n*                   JOBVS = 'V', VS contains the matrix which reduces A\n*                   to its partially converged Schur form.\n*             = N+1: the eigenvalues could not be reordered because some\n*                   eigenvalues were too close to separate (the problem\n*                   is very ill-conditioned);\n*             = N+2: after reordering, roundoff changed values of some\n*                   complex eigenvalues so that leading eigenvalues in\n*                   the Schur form no longer satisfy SELECT=.TRUE.  This\n*                   could also be caused by underflow due to scaling.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sdim, wr, wi, vs, work, info, a = NumRu::Lapack.dgees( jobvs, sort, a, lwork, [:usage => usage, :help => help]){|a,b| ... }\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_jobvs = argv[0];
  rblapack_sort = argv[1];
  rblapack_a = argv[2];
  rblapack_lwork = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  jobvs = StringValueCStr(rblapack_jobvs)[0];
  lwork = NUM2INT(rblapack_lwork);
  sort = StringValueCStr(rblapack_sort)[0];
  ldvs = lsame_(&jobvs,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_wr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rblapack_wr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_wi = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rblapack_wi, doublereal*);
  {
    int shape[2];
    shape[0] = ldvs;
    shape[1] = n;
    rblapack_vs = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vs = NA_PTR_TYPE(rblapack_vs, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  dgees_(&jobvs, &sort, rblapack_select, &n, a, &lda, &sdim, wr, wi, vs, &ldvs, work, &lwork, bwork, &info);

  free(bwork);
  rblapack_sdim = INT2NUM(sdim);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(7, rblapack_sdim, rblapack_wr, rblapack_wi, rblapack_vs, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_dgees(VALUE mLapack){
  rb_define_module_function(mLapack, "dgees", rblapack_dgees, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
