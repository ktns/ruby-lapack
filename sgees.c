#include "rb_lapack.h"

static logical
rb_select(real *arg0, real *arg1){
  VALUE rb_arg0, rb_arg1;

  VALUE rb_ret;
  logical ret;

  rb_arg0 = rb_float_new((double)(*arg0));
  rb_arg1 = rb_float_new((double)(*arg1));

  rb_ret = rb_yield_values(2, rb_arg0, rb_arg1);

  ret = (rb_ret == Qtrue);
  return ret;
}

extern VOID sgees_(char *jobvs, char *sort, L_fp *select, integer *n, real *a, integer *lda, integer *sdim, real *wr, real *wi, real *vs, integer *ldvs, real *work, integer *lwork, logical *bwork, integer *info);

static VALUE
rb_sgees(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvs;
  char jobvs; 
  VALUE rb_sort;
  char sort; 
  VALUE rb_a;
  real *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_sdim;
  integer sdim; 
  VALUE rb_wr;
  real *wr; 
  VALUE rb_wi;
  real *wi; 
  VALUE rb_vs;
  real *vs; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  logical *bwork;

  integer lda;
  integer n;
  integer ldvs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sdim, wr, wi, vs, work, info, a = NumRu::Lapack.sgees( jobvs, sort, a, lwork){|a,b| ... }\n    or\n  NumRu::Lapack.sgees  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI, VS, LDVS, WORK, LWORK, BWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGEES computes for an N-by-N real nonsymmetric matrix A, the\n*  eigenvalues, the real Schur form T, and, optionally, the matrix of\n*  Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).\n*\n*  Optionally, it also orders the eigenvalues on the diagonal of the\n*  real Schur form so that selected eigenvalues are at the top left.\n*  The leading columns of Z then form an orthonormal basis for the\n*  invariant subspace corresponding to the selected eigenvalues.\n*\n*  A matrix is in real Schur form if it is upper quasi-triangular with\n*  1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the\n*  form\n*          [  a  b  ]\n*          [  c  a  ]\n*\n*  where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVS   (input) CHARACTER*1\n*          = 'N': Schur vectors are not computed;\n*          = 'V': Schur vectors are computed.\n*\n*  SORT    (input) CHARACTER*1\n*          Specifies whether or not to order the eigenvalues on the\n*          diagonal of the Schur form.\n*          = 'N': Eigenvalues are not ordered;\n*          = 'S': Eigenvalues are ordered (see SELECT).\n*\n*  SELECT  (external procedure) LOGICAL FUNCTION of two REAL arguments\n*          SELECT must be declared EXTERNAL in the calling subroutine.\n*          If SORT = 'S', SELECT is used to select eigenvalues to sort\n*          to the top left of the Schur form.\n*          If SORT = 'N', SELECT is not referenced.\n*          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if\n*          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex\n*          conjugate pair of eigenvalues is selected, then both complex\n*          eigenvalues are selected.\n*          Note that a selected complex eigenvalue may no longer\n*          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since\n*          ordering may change the value of complex eigenvalues\n*          (especially if the eigenvalue is ill-conditioned); in this\n*          case INFO is set to N+2 (see INFO below).\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the N-by-N matrix A.\n*          On exit, A has been overwritten by its real Schur form T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  SDIM    (output) INTEGER\n*          If SORT = 'N', SDIM = 0.\n*          If SORT = 'S', SDIM = number of eigenvalues (after sorting)\n*                         for which SELECT is true. (Complex conjugate\n*                         pairs for which SELECT is true for either\n*                         eigenvalue count as 2.)\n*\n*  WR      (output) REAL array, dimension (N)\n*  WI      (output) REAL array, dimension (N)\n*          WR and WI contain the real and imaginary parts,\n*          respectively, of the computed eigenvalues in the same order\n*          that they appear on the diagonal of the output Schur form T.\n*          Complex conjugate pairs of eigenvalues will appear\n*          consecutively with the eigenvalue having the positive\n*          imaginary part first.\n*\n*  VS      (output) REAL array, dimension (LDVS,N)\n*          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur\n*          vectors.\n*          If JOBVS = 'N', VS is not referenced.\n*\n*  LDVS    (input) INTEGER\n*          The leading dimension of the array VS.  LDVS >= 1; if\n*          JOBVS = 'V', LDVS >= N.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) contains the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,3*N).\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  BWORK   (workspace) LOGICAL array, dimension (N)\n*          Not referenced if SORT = 'N'.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value.\n*          > 0: if INFO = i, and i is\n*             <= N: the QR algorithm failed to compute all the\n*                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI\n*                   contain those eigenvalues which have converged; if\n*                   JOBVS = 'V', VS contains the matrix which reduces A\n*                   to its partially converged Schur form.\n*             = N+1: the eigenvalues could not be reordered because some\n*                   eigenvalues were too close to separate (the problem\n*                   is very ill-conditioned);\n*             = N+2: after reordering, roundoff changed values of some\n*                   complex eigenvalues so that leading eigenvalues in\n*                   the Schur form no longer satisfy SELECT=.TRUE.  This\n*                   could also be caused by underflow due to scaling.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_jobvs = argv[0];
  rb_sort = argv[1];
  rb_a = argv[2];
  rb_lwork = argv[3];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  jobvs = StringValueCStr(rb_jobvs)[0];
  lwork = NUM2INT(rb_lwork);
  sort = StringValueCStr(rb_sort)[0];
  ldvs = lsame_(&jobvs,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_wr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_wi = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, real*);
  {
    int shape[2];
    shape[0] = ldvs;
    shape[1] = n;
    rb_vs = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vs = NA_PTR_TYPE(rb_vs, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
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
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  sgees_(&jobvs, &sort, rb_select, &n, a, &lda, &sdim, wr, wi, vs, &ldvs, work, &lwork, bwork, &info);

  free(bwork);
  rb_sdim = INT2NUM(sdim);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_sdim, rb_wr, rb_wi, rb_vs, rb_work, rb_info, rb_a);
}

void
init_lapack_sgees(VALUE mLapack){
  rb_define_module_function(mLapack, "sgees", rb_sgees, -1);
}
