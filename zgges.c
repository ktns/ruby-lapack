#include "rb_lapack.h"

static logical
rb_selctg(doublecomplex *arg0, doublecomplex *arg1){
  VALUE rb_arg0, rb_arg1;

  VALUE rb_ret;
  logical ret;

  rb_arg0 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(arg0->r)), rb_float_new((double)(arg0->i)));
  rb_arg1 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(arg1->r)), rb_float_new((double)(arg1->i)));

  rb_ret = rb_yield_values(2, rb_arg0, rb_arg1);

  ret = (rb_ret == Qtrue);
  return ret;
}

extern VOID zgges_(char *jobvsl, char *jobvsr, char *sort, L_fp *selctg, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, integer *sdim, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, integer *ldvsl, doublecomplex *vsr, integer *ldvsr, doublecomplex *work, integer *lwork, doublereal *rwork, logical *bwork, integer *info);

static VALUE
rb_zgges(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvsl;
  char jobvsl; 
  VALUE rb_jobvsr;
  char jobvsr; 
  VALUE rb_sort;
  char sort; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_sdim;
  integer sdim; 
  VALUE rb_alpha;
  doublecomplex *alpha; 
  VALUE rb_beta;
  doublecomplex *beta; 
  VALUE rb_vsl;
  doublecomplex *vsl; 
  VALUE rb_vsr;
  doublecomplex *vsr; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  doublereal *rwork;
  logical *bwork;

  integer lda;
  integer n;
  integer ldb;
  integer ldvsl;
  integer ldvsr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sdim, alpha, beta, vsl, vsr, work, info, a, b = NumRu::Lapack.zgges( jobvsl, jobvsr, sort, a, b, lwork){|a,b| ... }\n    or\n  NumRu::Lapack.zgges  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, BWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGGES computes for a pair of N-by-N complex nonsymmetric matrices\n*  (A,B), the generalized eigenvalues, the generalized complex Schur\n*  form (S, T), and optionally left and/or right Schur vectors (VSL\n*  and VSR). This gives the generalized Schur factorization\n*\n*          (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )\n*\n*  where (VSR)**H is the conjugate-transpose of VSR.\n*\n*  Optionally, it also orders the eigenvalues so that a selected cluster\n*  of eigenvalues appears in the leading diagonal blocks of the upper\n*  triangular matrix S and the upper triangular matrix T. The leading\n*  columns of VSL and VSR then form an unitary basis for the\n*  corresponding left and right eigenspaces (deflating subspaces).\n*\n*  (If only the generalized eigenvalues are needed, use the driver\n*  ZGGEV instead, which is faster.)\n*\n*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar w\n*  or a ratio alpha/beta = w, such that  A - w*B is singular.  It is\n*  usually represented as the pair (alpha,beta), as there is a\n*  reasonable interpretation for beta=0, and even for both being zero.\n*\n*  A pair of matrices (S,T) is in generalized complex Schur form if S\n*  and T are upper triangular and, in addition, the diagonal elements\n*  of T are non-negative real numbers.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVSL  (input) CHARACTER*1\n*          = 'N':  do not compute the left Schur vectors;\n*          = 'V':  compute the left Schur vectors.\n*\n*  JOBVSR  (input) CHARACTER*1\n*          = 'N':  do not compute the right Schur vectors;\n*          = 'V':  compute the right Schur vectors.\n*\n*  SORT    (input) CHARACTER*1\n*          Specifies whether or not to order the eigenvalues on the\n*          diagonal of the generalized Schur form.\n*          = 'N':  Eigenvalues are not ordered;\n*          = 'S':  Eigenvalues are ordered (see SELCTG).\n*\n*  SELCTG  (external procedure) LOGICAL FUNCTION of two COMPLEX*16 arguments\n*          SELCTG must be declared EXTERNAL in the calling subroutine.\n*          If SORT = 'N', SELCTG is not referenced.\n*          If SORT = 'S', SELCTG is used to select eigenvalues to sort\n*          to the top left of the Schur form.\n*          An eigenvalue ALPHA(j)/BETA(j) is selected if\n*          SELCTG(ALPHA(j),BETA(j)) is true.\n*\n*          Note that a selected complex eigenvalue may no longer satisfy\n*          SELCTG(ALPHA(j),BETA(j)) = .TRUE. after ordering, since\n*          ordering may change the value of complex eigenvalues\n*          (especially if the eigenvalue is ill-conditioned), in this\n*          case INFO is set to N+2 (See INFO below).\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VSL, and VSR.  N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)\n*          On entry, the first of the pair of matrices.\n*          On exit, A has been overwritten by its generalized Schur\n*          form S.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB, N)\n*          On entry, the second of the pair of matrices.\n*          On exit, B has been overwritten by its generalized Schur\n*          form T.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  SDIM    (output) INTEGER\n*          If SORT = 'N', SDIM = 0.\n*          If SORT = 'S', SDIM = number of eigenvalues (after sorting)\n*          for which SELCTG is true.\n*\n*  ALPHA   (output) COMPLEX*16 array, dimension (N)\n*  BETA    (output) COMPLEX*16 array, dimension (N)\n*          On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the\n*          generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j),\n*          j=1,...,N  are the diagonals of the complex Schur form (A,B)\n*          output by ZGGES. The  BETA(j) will be non-negative real.\n*\n*          Note: the quotients ALPHA(j)/BETA(j) may easily over- or\n*          underflow, and BETA(j) may even be zero.  Thus, the user\n*          should avoid naively computing the ratio alpha/beta.\n*          However, ALPHA will be always less than and usually\n*          comparable with norm(A) in magnitude, and BETA always less\n*          than and usually comparable with norm(B).\n*\n*  VSL     (output) COMPLEX*16 array, dimension (LDVSL,N)\n*          If JOBVSL = 'V', VSL will contain the left Schur vectors.\n*          Not referenced if JOBVSL = 'N'.\n*\n*  LDVSL   (input) INTEGER\n*          The leading dimension of the matrix VSL. LDVSL >= 1, and\n*          if JOBVSL = 'V', LDVSL >= N.\n*\n*  VSR     (output) COMPLEX*16 array, dimension (LDVSR,N)\n*          If JOBVSR = 'V', VSR will contain the right Schur vectors.\n*          Not referenced if JOBVSR = 'N'.\n*\n*  LDVSR   (input) INTEGER\n*          The leading dimension of the matrix VSR. LDVSR >= 1, and\n*          if JOBVSR = 'V', LDVSR >= N.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,2*N).\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (8*N)\n*\n*  BWORK   (workspace) LOGICAL array, dimension (N)\n*          Not referenced if SORT = 'N'.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          =1,...,N:\n*                The QZ iteration failed.  (A,B) are not in Schur\n*                form, but ALPHA(j) and BETA(j) should be correct for\n*                j=INFO+1,...,N.\n*          > N:  =N+1: other than QZ iteration failed in ZHGEQZ\n*                =N+2: after reordering, roundoff changed values of\n*                      some complex eigenvalues so that leading\n*                      eigenvalues in the Generalized Schur form no\n*                      longer satisfy SELCTG=.TRUE.  This could also\n*                      be caused due to scaling.\n*                =N+3: reordering falied in ZTGSEN.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_jobvsl = argv[0];
  rb_jobvsr = argv[1];
  rb_sort = argv[2];
  rb_a = argv[3];
  rb_b = argv[4];
  rb_lwork = argv[5];

  jobvsl = StringValueCStr(rb_jobvsl)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  jobvsr = StringValueCStr(rb_jobvsr)[0];
  sort = StringValueCStr(rb_sort)[0];
  lwork = NUM2INT(rb_lwork);
  ldvsr = lsame_(&jobvsr,"V") ? n : 1;
  ldvsl = lsame_(&jobvsl,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_alpha = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  alpha = NA_PTR_TYPE(rb_alpha, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_beta = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldvsl;
    shape[1] = n;
    rb_vsl = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vsl = NA_PTR_TYPE(rb_vsl, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldvsr;
    shape[1] = n;
    rb_vsr = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vsr = NA_PTR_TYPE(rb_vsr, doublecomplex*);
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
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  rwork = ALLOC_N(doublereal, (8*n));
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  zgges_(&jobvsl, &jobvsr, &sort, rb_selctg, &n, a, &lda, b, &ldb, &sdim, alpha, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, rwork, bwork, &info);

  free(rwork);
  free(bwork);
  rb_sdim = INT2NUM(sdim);
  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_sdim, rb_alpha, rb_beta, rb_vsl, rb_vsr, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_zgges(VALUE mLapack){
  rb_define_module_function(mLapack, "zgges", rb_zgges, -1);
}
