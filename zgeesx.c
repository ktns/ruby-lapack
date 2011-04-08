#include "rb_lapack.h"

static logical
rblapack_select(doublecomplex *arg0){
  VALUE rblapack_arg0;

  VALUE rblapack_ret;
  logical ret;

  rblapack_arg0 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(arg0->r)), rb_float_new((double)(arg0->i)));

  rblapack_ret = rb_yield_values(1, rblapack_arg0);

  ret = (rblapack_ret == Qtrue);
  return ret;
}

extern VOID zgeesx_(char *jobvs, char *sort, L_fp *select, char *sense, integer *n, doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w, doublecomplex *vs, integer *ldvs, doublereal *rconde, doublereal *rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, logical *bwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgeesx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobvs;
  char jobvs; 
  VALUE rblapack_sort;
  char sort; 
  VALUE rblapack_sense;
  char sense; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_sdim;
  integer sdim; 
  VALUE rblapack_w;
  doublecomplex *w; 
  VALUE rblapack_vs;
  doublecomplex *vs; 
  VALUE rblapack_rconde;
  doublereal rconde; 
  VALUE rblapack_rcondv;
  doublereal rcondv; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;
  doublereal *rwork;
  logical *bwork;

  integer lda;
  integer n;
  integer ldvs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sdim, w, vs, rconde, rcondv, work, info, a = NumRu::Lapack.zgeesx( jobvs, sort, sense, a, lwork, [:usage => usage, :help => help]){|a| ... }\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, W, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGEESX computes for an N-by-N complex nonsymmetric matrix A, the\n*  eigenvalues, the Schur form T, and, optionally, the matrix of Schur\n*  vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).\n*\n*  Optionally, it also orders the eigenvalues on the diagonal of the\n*  Schur form so that selected eigenvalues are at the top left;\n*  computes a reciprocal condition number for the average of the\n*  selected eigenvalues (RCONDE); and computes a reciprocal condition\n*  number for the right invariant subspace corresponding to the\n*  selected eigenvalues (RCONDV).  The leading columns of Z form an\n*  orthonormal basis for this invariant subspace.\n*\n*  For further explanation of the reciprocal condition numbers RCONDE\n*  and RCONDV, see Section 4.10 of the LAPACK Users' Guide (where\n*  these quantities are called s and sep respectively).\n*\n*  A complex matrix is in Schur form if it is upper triangular.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVS   (input) CHARACTER*1\n*          = 'N': Schur vectors are not computed;\n*          = 'V': Schur vectors are computed.\n*\n*  SORT    (input) CHARACTER*1\n*          Specifies whether or not to order the eigenvalues on the\n*          diagonal of the Schur form.\n*          = 'N': Eigenvalues are not ordered;\n*          = 'S': Eigenvalues are ordered (see SELECT).\n*\n*  SELECT  (external procedure) LOGICAL FUNCTION of one COMPLEX*16 argument\n*          SELECT must be declared EXTERNAL in the calling subroutine.\n*          If SORT = 'S', SELECT is used to select eigenvalues to order\n*          to the top left of the Schur form.\n*          If SORT = 'N', SELECT is not referenced.\n*          An eigenvalue W(j) is selected if SELECT(W(j)) is true.\n*\n*  SENSE   (input) CHARACTER*1\n*          Determines which reciprocal condition numbers are computed.\n*          = 'N': None are computed;\n*          = 'E': Computed for average of selected eigenvalues only;\n*          = 'V': Computed for selected right invariant subspace only;\n*          = 'B': Computed for both.\n*          If SENSE = 'E', 'V' or 'B', SORT must equal 'S'.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)\n*          On entry, the N-by-N matrix A.\n*          On exit, A is overwritten by its Schur form T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  SDIM    (output) INTEGER\n*          If SORT = 'N', SDIM = 0.\n*          If SORT = 'S', SDIM = number of eigenvalues for which\n*                         SELECT is true.\n*\n*  W       (output) COMPLEX*16 array, dimension (N)\n*          W contains the computed eigenvalues, in the same order\n*          that they appear on the diagonal of the output Schur form T.\n*\n*  VS      (output) COMPLEX*16 array, dimension (LDVS,N)\n*          If JOBVS = 'V', VS contains the unitary matrix Z of Schur\n*          vectors.\n*          If JOBVS = 'N', VS is not referenced.\n*\n*  LDVS    (input) INTEGER\n*          The leading dimension of the array VS.  LDVS >= 1, and if\n*          JOBVS = 'V', LDVS >= N.\n*\n*  RCONDE  (output) DOUBLE PRECISION\n*          If SENSE = 'E' or 'B', RCONDE contains the reciprocal\n*          condition number for the average of the selected eigenvalues.\n*          Not referenced if SENSE = 'N' or 'V'.\n*\n*  RCONDV  (output) DOUBLE PRECISION\n*          If SENSE = 'V' or 'B', RCONDV contains the reciprocal\n*          condition number for the selected right invariant subspace.\n*          Not referenced if SENSE = 'N' or 'E'.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,2*N).\n*          Also, if SENSE = 'E' or 'V' or 'B', LWORK >= 2*SDIM*(N-SDIM),\n*          where SDIM is the number of selected eigenvalues computed by\n*          this routine.  Note that 2*SDIM*(N-SDIM) <= N*N/2. Note also\n*          that an error is only returned if LWORK < max(1,2*N), but if\n*          SENSE = 'E' or 'V' or 'B' this may not be large enough.\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates upper bound on the optimal size of the\n*          array WORK, returns this value as the first entry of the WORK\n*          array, and no error message related to LWORK is issued by\n*          XERBLA.\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  BWORK   (workspace) LOGICAL array, dimension (N)\n*          Not referenced if SORT = 'N'.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value.\n*          > 0: if INFO = i, and i is\n*             <= N: the QR algorithm failed to compute all the\n*                   eigenvalues; elements 1:ILO-1 and i+1:N of W\n*                   contain those eigenvalues which have converged; if\n*                   JOBVS = 'V', VS contains the transformation which\n*                   reduces A to its partially converged Schur form.\n*             = N+1: the eigenvalues could not be reordered because some\n*                   eigenvalues were too close to separate (the problem\n*                   is very ill-conditioned);\n*             = N+2: after reordering, roundoff changed values of some\n*                   complex eigenvalues so that leading eigenvalues in\n*                   the Schur form no longer satisfy SELECT=.TRUE.  This\n*                   could also be caused by underflow due to scaling.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sdim, w, vs, rconde, rcondv, work, info, a = NumRu::Lapack.zgeesx( jobvs, sort, sense, a, lwork, [:usage => usage, :help => help]){|a| ... }\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_jobvs = argv[0];
  rblapack_sort = argv[1];
  rblapack_sense = argv[2];
  rblapack_a = argv[3];
  rblapack_lwork = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  jobvs = StringValueCStr(rblapack_jobvs)[0];
  sort = StringValueCStr(rblapack_sort)[0];
  lwork = NUM2INT(rblapack_lwork);
  sense = StringValueCStr(rblapack_sense)[0];
  ldvs = lsame_(&jobvs,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldvs;
    shape[1] = n;
    rblapack_vs = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vs = NA_PTR_TYPE(rblapack_vs, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  rwork = ALLOC_N(doublereal, (n));
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  zgeesx_(&jobvs, &sort, rblapack_select, &sense, &n, a, &lda, &sdim, w, vs, &ldvs, &rconde, &rcondv, work, &lwork, rwork, bwork, &info);

  free(rwork);
  free(bwork);
  rblapack_sdim = INT2NUM(sdim);
  rblapack_rconde = rb_float_new((double)rconde);
  rblapack_rcondv = rb_float_new((double)rcondv);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(8, rblapack_sdim, rblapack_w, rblapack_vs, rblapack_rconde, rblapack_rcondv, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_zgeesx(VALUE mLapack){
  rb_define_module_function(mLapack, "zgeesx", rblapack_zgeesx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
