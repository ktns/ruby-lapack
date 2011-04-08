#include "rb_lapack.h"

static logical
rblapack_select(complex *arg0){
  VALUE rblapack_arg0;

  VALUE rblapack_ret;
  logical ret;

  rblapack_arg0 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(arg0->r)), rb_float_new((double)(arg0->i)));

  rblapack_ret = rb_yield_values(1, rblapack_arg0);

  ret = (rblapack_ret == Qtrue);
  return ret;
}

extern VOID cgees_(char *jobvs, char *sort, L_fp *select, integer *n, complex *a, integer *lda, integer *sdim, complex *w, complex *vs, integer *ldvs, complex *work, integer *lwork, real *rwork, logical *bwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgees(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobvs;
  char jobvs; 
  VALUE rblapack_sort;
  char sort; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_sdim;
  integer sdim; 
  VALUE rblapack_w;
  complex *w; 
  VALUE rblapack_vs;
  complex *vs; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  real *rwork;
  logical *bwork;

  integer lda;
  integer n;
  integer ldvs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sdim, w, vs, work, info, a = NumRu::Lapack.cgees( jobvs, sort, a, lwork, [:usage => usage, :help => help]){|a| ... }\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS, LDVS, WORK, LWORK, RWORK, BWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGEES computes for an N-by-N complex nonsymmetric matrix A, the\n*  eigenvalues, the Schur form T, and, optionally, the matrix of Schur\n*  vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).\n*\n*  Optionally, it also orders the eigenvalues on the diagonal of the\n*  Schur form so that selected eigenvalues are at the top left.\n*  The leading columns of Z then form an orthonormal basis for the\n*  invariant subspace corresponding to the selected eigenvalues.\n\n*  A complex matrix is in Schur form if it is upper triangular.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVS   (input) CHARACTER*1\n*          = 'N': Schur vectors are not computed;\n*          = 'V': Schur vectors are computed.\n*\n*  SORT    (input) CHARACTER*1\n*          Specifies whether or not to order the eigenvalues on the\n*          diagonal of the Schur form.\n*          = 'N': Eigenvalues are not ordered:\n*          = 'S': Eigenvalues are ordered (see SELECT).\n*\n*  SELECT  (external procedure) LOGICAL FUNCTION of one COMPLEX argument\n*          SELECT must be declared EXTERNAL in the calling subroutine.\n*          If SORT = 'S', SELECT is used to select eigenvalues to order\n*          to the top left of the Schur form.\n*          IF SORT = 'N', SELECT is not referenced.\n*          The eigenvalue W(j) is selected if SELECT(W(j)) is true.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the N-by-N matrix A.\n*          On exit, A has been overwritten by its Schur form T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  SDIM    (output) INTEGER\n*          If SORT = 'N', SDIM = 0.\n*          If SORT = 'S', SDIM = number of eigenvalues for which\n*                         SELECT is true.\n*\n*  W       (output) COMPLEX array, dimension (N)\n*          W contains the computed eigenvalues, in the same order that\n*          they appear on the diagonal of the output Schur form T.\n*\n*  VS      (output) COMPLEX array, dimension (LDVS,N)\n*          If JOBVS = 'V', VS contains the unitary matrix Z of Schur\n*          vectors.\n*          If JOBVS = 'N', VS is not referenced.\n*\n*  LDVS    (input) INTEGER\n*          The leading dimension of the array VS.  LDVS >= 1; if\n*          JOBVS = 'V', LDVS >= N.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,2*N).\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) REAL array, dimension (N)\n*\n*  BWORK   (workspace) LOGICAL array, dimension (N)\n*          Not referenced if SORT = 'N'.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value.\n*          > 0: if INFO = i, and i is\n*               <= N:  the QR algorithm failed to compute all the\n*                      eigenvalues; elements 1:ILO-1 and i+1:N of W\n*                      contain those eigenvalues which have converged;\n*                      if JOBVS = 'V', VS contains the matrix which\n*                      reduces A to its partially converged Schur form.\n*               = N+1: the eigenvalues could not be reordered because\n*                      some eigenvalues were too close to separate (the\n*                      problem is very ill-conditioned);\n*               = N+2: after reordering, roundoff changed values of\n*                      some complex eigenvalues so that leading\n*                      eigenvalues in the Schur form no longer satisfy\n*                      SELECT = .TRUE..  This could also be caused by\n*                      underflow due to scaling.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sdim, w, vs, work, info, a = NumRu::Lapack.cgees( jobvs, sort, a, lwork, [:usage => usage, :help => help]){|a| ... }\n");
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
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  jobvs = StringValueCStr(rblapack_jobvs)[0];
  lwork = NUM2INT(rblapack_lwork);
  sort = StringValueCStr(rblapack_sort)[0];
  ldvs = lsame_(&jobvs,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, complex*);
  {
    int shape[2];
    shape[0] = ldvs;
    shape[1] = n;
    rblapack_vs = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vs = NA_PTR_TYPE(rblapack_vs, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  rwork = ALLOC_N(real, (n));
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  cgees_(&jobvs, &sort, rblapack_select, &n, a, &lda, &sdim, w, vs, &ldvs, work, &lwork, rwork, bwork, &info);

  free(rwork);
  free(bwork);
  rblapack_sdim = INT2NUM(sdim);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_sdim, rblapack_w, rblapack_vs, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_cgees(VALUE mLapack){
  rb_define_module_function(mLapack, "cgees", rblapack_cgees, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
