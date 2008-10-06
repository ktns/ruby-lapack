#include "rb_lapack.h"

static logical
rb_selctg(real *arg0, real *arg1, real *arg2){
  VALUE rb_arg0, rb_arg1, rb_arg2;

  VALUE rb_ret;
  logical ret;

  rb_arg0 = rb_float_new((double)(*arg0));
  rb_arg1 = rb_float_new((double)(*arg1));
  rb_arg2 = rb_float_new((double)(*arg2));

  rb_ret = rb_yield_values(3, rb_arg0, rb_arg1, rb_arg2);

  ret = (rb_ret == Qtrue);
  return ret;
}

static VALUE
rb_sgges(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvsl;
  char jobvsl; 
  VALUE rb_jobvsr;
  char jobvsr; 
  VALUE rb_sort;
  char sort; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_sdim;
  integer sdim; 
  VALUE rb_alphar;
  real *alphar; 
  VALUE rb_alphai;
  real *alphai; 
  VALUE rb_beta;
  real *beta; 
  VALUE rb_vsl;
  real *vsl; 
  VALUE rb_vsr;
  real *vsr; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  logical *bwork;

  integer lda;
  integer n;
  integer ldb;
  integer ldvsl;
  integer ldvsr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sdim, alphar, alphai, beta, vsl, vsr, work, info, a, b = NumRu::Lapack.sgges( jobvsl, jobvsr, sort, a, b, lwork){|a,b,c| ... }\n    or\n  NumRu::Lapack.sgges  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, BWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGGES computes for a pair of N-by-N real nonsymmetric matrices (A,B),\n*  the generalized eigenvalues, the generalized real Schur form (S,T),\n*  optionally, the left and/or right matrices of Schur vectors (VSL and\n*  VSR). This gives the generalized Schur factorization\n*\n*           (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )\n*\n*  Optionally, it also orders the eigenvalues so that a selected cluster\n*  of eigenvalues appears in the leading diagonal blocks of the upper\n*  quasi-triangular matrix S and the upper triangular matrix T.The\n*  leading columns of VSL and VSR then form an orthonormal basis for the\n*  corresponding left and right eigenspaces (deflating subspaces).\n*\n*  (If only the generalized eigenvalues are needed, use the driver\n*  SGGEV instead, which is faster.)\n*\n*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar w\n*  or a ratio alpha/beta = w, such that  A - w*B is singular.  It is\n*  usually represented as the pair (alpha,beta), as there is a\n*  reasonable interpretation for beta=0 or both being zero.\n*\n*  A pair of matrices (S,T) is in generalized real Schur form if T is\n*  upper triangular with non-negative diagonal and S is block upper\n*  triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond\n*  to real generalized eigenvalues, while 2-by-2 blocks of S will be\n*  \"standardized\" by making the corresponding elements of T have the\n*  form:\n*          [  a  0  ]\n*          [  0  b  ]\n*\n*  and the pair of corresponding 2-by-2 blocks in S and T will have a\n*  complex conjugate pair of generalized eigenvalues.\n*\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVSL  (input) CHARACTER*1\n*          = 'N':  do not compute the left Schur vectors;\n*          = 'V':  compute the left Schur vectors.\n*\n*  JOBVSR  (input) CHARACTER*1\n*          = 'N':  do not compute the right Schur vectors;\n*          = 'V':  compute the right Schur vectors.\n*\n*  SORT    (input) CHARACTER*1\n*          Specifies whether or not to order the eigenvalues on the\n*          diagonal of the generalized Schur form.\n*          = 'N':  Eigenvalues are not ordered;\n*          = 'S':  Eigenvalues are ordered (see SELCTG);\n*\n*  SELCTG  (external procedure) LOGICAL FUNCTION of three REAL arguments\n*          SELCTG must be declared EXTERNAL in the calling subroutine.\n*          If SORT = 'N', SELCTG is not referenced.\n*          If SORT = 'S', SELCTG is used to select eigenvalues to sort\n*          to the top left of the Schur form.\n*          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if\n*          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either\n*          one of a complex conjugate pair of eigenvalues is selected,\n*          then both complex eigenvalues are selected.\n*\n*          Note that in the ill-conditioned case, a selected complex\n*          eigenvalue may no longer satisfy SELCTG(ALPHAR(j),ALPHAI(j),\n*          BETA(j)) = .TRUE. after ordering. INFO is to be set to N+2\n*          in this case.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VSL, and VSR.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA, N)\n*          On entry, the first of the pair of matrices.\n*          On exit, A has been overwritten by its generalized Schur\n*          form S.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) REAL array, dimension (LDB, N)\n*          On entry, the second of the pair of matrices.\n*          On exit, B has been overwritten by its generalized Schur\n*          form T.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  SDIM    (output) INTEGER\n*          If SORT = 'N', SDIM = 0.\n*          If SORT = 'S', SDIM = number of eigenvalues (after sorting)\n*          for which SELCTG is true.  (Complex conjugate pairs for which\n*          SELCTG is true for either eigenvalue count as 2.)\n*\n*  ALPHAR  (output) REAL array, dimension (N)\n*  ALPHAI  (output) REAL array, dimension (N)\n*  BETA    (output) REAL array, dimension (N)\n*          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will\n*          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,\n*          and  BETA(j),j=1,...,N are the diagonals of the complex Schur\n*          form (S,T) that would result if the 2-by-2 diagonal blocks of\n*          the real Schur form of (A,B) were further reduced to\n*          triangular form using 2-by-2 complex unitary transformations.\n*          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if\n*          positive, then the j-th and (j+1)-st eigenvalues are a\n*          complex conjugate pair, with ALPHAI(j+1) negative.\n*\n*          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)\n*          may easily over- or underflow, and BETA(j) may even be zero.\n*          Thus, the user should avoid naively computing the ratio.\n*          However, ALPHAR and ALPHAI will be always less than and\n*          usually comparable with norm(A) in magnitude, and BETA always\n*          less than and usually comparable with norm(B).\n*\n*  VSL     (output) REAL array, dimension (LDVSL,N)\n*          If JOBVSL = 'V', VSL will contain the left Schur vectors.\n*          Not referenced if JOBVSL = 'N'.\n*\n*  LDVSL   (input) INTEGER\n*          The leading dimension of the matrix VSL. LDVSL >=1, and\n*          if JOBVSL = 'V', LDVSL >= N.\n*\n*  VSR     (output) REAL array, dimension (LDVSR,N)\n*          If JOBVSR = 'V', VSR will contain the right Schur vectors.\n*          Not referenced if JOBVSR = 'N'.\n*\n*  LDVSR   (input) INTEGER\n*          The leading dimension of the matrix VSR. LDVSR >= 1, and\n*          if JOBVSR = 'V', LDVSR >= N.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If N = 0, LWORK >= 1, else LWORK >= max(8*N,6*N+16).\n*          For good performance , LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  BWORK   (workspace) LOGICAL array, dimension (N)\n*          Not referenced if SORT = 'N'.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1,...,N:\n*                The QZ iteration failed.  (A,B) are not in Schur\n*                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should\n*                be correct for j=INFO+1,...,N.\n*          > N:  =N+1: other than QZ iteration failed in SHGEQZ.\n*                =N+2: after reordering, roundoff changed values of\n*                      some complex eigenvalues so that leading\n*                      eigenvalues in the Generalized Schur form no\n*                      longer satisfy SELCTG=.TRUE.  This could also\n*                      be caused due to scaling.\n*                =N+3: reordering failed in STGSEN.\n*\n\n*  =====================================================================\n*\n\n");
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
  jobvsr = StringValueCStr(rb_jobvsr)[0];
  sort = StringValueCStr(rb_sort)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_alphar = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alphar = NA_PTR_TYPE(rb_alphar, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_alphai = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alphai = NA_PTR_TYPE(rb_alphai, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_beta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, real*);
  ldvsl = lsame_(&jobvsl,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldvsl;
    shape[1] = n;
    rb_vsl = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vsl = NA_PTR_TYPE(rb_vsl, real*);
  ldvsr = lsame_(&jobvsr,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldvsr;
    shape[1] = n;
    rb_vsr = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vsr = NA_PTR_TYPE(rb_vsr, real*);
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
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  sgges_(&jobvsl, &jobvsr, &sort, rb_selctg, &n, a, &lda, b, &ldb, &sdim, alphar, alphai, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, bwork, &info);

  free(bwork);
  rb_sdim = INT2NUM(sdim);
  rb_info = INT2NUM(info);
  return rb_ary_new3(10, rb_sdim, rb_alphar, rb_alphai, rb_beta, rb_vsl, rb_vsr, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_sgges(VALUE mLapack){
  rb_define_module_function(mLapack, "sgges", rb_sgges, -1);
}
