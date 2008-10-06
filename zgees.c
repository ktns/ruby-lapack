#include "rb_lapack.h"

static logical
rb_select(doublecomplex *arg0){
  VALUE rb_arg0;

  VALUE rb_ret;
  logical ret;

  rb_arg0 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(arg0->r)), rb_float_new((double)(arg0->i)));

  rb_ret = rb_yield_values(1, rb_arg0);

  ret = (rb_ret == Qtrue);
  return ret;
}

static VALUE
rb_zgees(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvs;
  char jobvs; 
  VALUE rb_sort;
  char sort; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_sdim;
  integer sdim; 
  VALUE rb_w;
  doublecomplex *w; 
  VALUE rb_vs;
  doublecomplex *vs; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  doublereal *rwork;
  logical *bwork;

  integer lda;
  integer n;
  integer ldvs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sdim, w, vs, work, info, a = NumRu::Lapack.zgees( jobvs, sort, a, lwork){|a| ... }\n    or\n  NumRu::Lapack.zgees  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS, LDVS, WORK, LWORK, RWORK, BWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGEES computes for an N-by-N complex nonsymmetric matrix A, the\n*  eigenvalues, the Schur form T, and, optionally, the matrix of Schur\n*  vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).\n*\n*  Optionally, it also orders the eigenvalues on the diagonal of the\n*  Schur form so that selected eigenvalues are at the top left.\n*  The leading columns of Z then form an orthonormal basis for the\n*  invariant subspace corresponding to the selected eigenvalues.\n*\n*  A complex matrix is in Schur form if it is upper triangular.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVS   (input) CHARACTER*1\n*          = 'N': Schur vectors are not computed;\n*          = 'V': Schur vectors are computed.\n*\n*  SORT    (input) CHARACTER*1\n*          Specifies whether or not to order the eigenvalues on the\n*          diagonal of the Schur form.\n*          = 'N': Eigenvalues are not ordered:\n*          = 'S': Eigenvalues are ordered (see SELECT).\n*\n*  SELECT  (external procedure) LOGICAL FUNCTION of one COMPLEX*16 argument\n*          SELECT must be declared EXTERNAL in the calling subroutine.\n*          If SORT = 'S', SELECT is used to select eigenvalues to order\n*          to the top left of the Schur form.\n*          IF SORT = 'N', SELECT is not referenced.\n*          The eigenvalue W(j) is selected if SELECT(W(j)) is true.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the N-by-N matrix A.\n*          On exit, A has been overwritten by its Schur form T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  SDIM    (output) INTEGER\n*          If SORT = 'N', SDIM = 0.\n*          If SORT = 'S', SDIM = number of eigenvalues for which\n*                         SELECT is true.\n*\n*  W       (output) COMPLEX*16 array, dimension (N)\n*          W contains the computed eigenvalues, in the same order that\n*          they appear on the diagonal of the output Schur form T.\n*\n*  VS      (output) COMPLEX*16 array, dimension (LDVS,N)\n*          If JOBVS = 'V', VS contains the unitary matrix Z of Schur\n*          vectors.\n*          If JOBVS = 'N', VS is not referenced.\n*\n*  LDVS    (input) INTEGER\n*          The leading dimension of the array VS.  LDVS >= 1; if\n*          JOBVS = 'V', LDVS >= N.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,2*N).\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  BWORK   (workspace) LOGICAL array, dimension (N)\n*          Not referenced if SORT = 'N'.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value.\n*          > 0: if INFO = i, and i is\n*               <= N:  the QR algorithm failed to compute all the\n*                      eigenvalues; elements 1:ILO-1 and i+1:N of W\n*                      contain those eigenvalues which have converged;\n*                      if JOBVS = 'V', VS contains the matrix which\n*                      reduces A to its partially converged Schur form.\n*               = N+1: the eigenvalues could not be reordered because\n*                      some eigenvalues were too close to separate (the\n*                      problem is very ill-conditioned);\n*               = N+2: after reordering, roundoff changed values of\n*                      some complex eigenvalues so that leading\n*                      eigenvalues in the Schur form no longer satisfy\n*                      SELECT = .TRUE..  This could also be caused by\n*                      underflow due to scaling.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_jobvs = argv[0];
  rb_sort = argv[1];
  rb_a = argv[2];
  rb_lwork = argv[3];

  jobvs = StringValueCStr(rb_jobvs)[0];
  sort = StringValueCStr(rb_sort)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublecomplex*);
  ldvs = lsame_(&jobvs,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldvs;
    shape[1] = n;
    rb_vs = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vs = NA_PTR_TYPE(rb_vs, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  rwork = ALLOC_N(doublereal, (n));
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  zgees_(&jobvs, &sort, rb_select, &n, a, &lda, &sdim, w, vs, &ldvs, work, &lwork, rwork, bwork, &info);

  free(rwork);
  free(bwork);
  rb_sdim = INT2NUM(sdim);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_sdim, rb_w, rb_vs, rb_work, rb_info, rb_a);
}

void
init_lapack_zgees(VALUE mLapack){
  rb_define_module_function(mLapack, "zgees", rb_zgees, -1);
}
