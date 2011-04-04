#include "rb_lapack.h"

extern VOID sgegs_(char *jobvsl, char *jobvsr, integer *n, real *a, integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real *beta, real *vsl, integer *ldvsl, real *vsr, integer *ldvsr, real *work, integer *lwork, integer *info);

static VALUE
rb_sgegs(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvsl;
  char jobvsl; 
  VALUE rb_jobvsr;
  char jobvsr; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_lwork;
  integer lwork; 
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

  integer lda;
  integer n;
  integer ldb;
  integer ldvsl;
  integer ldvsr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alphar, alphai, beta, vsl, vsr, work, info, a, b = NumRu::Lapack.sgegs( jobvsl, jobvsr, a, b, lwork)\n    or\n  NumRu::Lapack.sgegs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine SGGES.\n*\n*  SGEGS computes the eigenvalues, real Schur form, and, optionally,\n*  left and or/right Schur vectors of a real matrix pair (A,B).\n*  Given two square matrices A and B, the generalized real Schur\n*  factorization has the form\n*  \n*    A = Q*S*Z**T,  B = Q*T*Z**T\n*\n*  where Q and Z are orthogonal matrices, T is upper triangular, and S\n*  is an upper quasi-triangular matrix with 1-by-1 and 2-by-2 diagonal\n*  blocks, the 2-by-2 blocks corresponding to complex conjugate pairs\n*  of eigenvalues of (A,B).  The columns of Q are the left Schur vectors\n*  and the columns of Z are the right Schur vectors.\n*  \n*  If only the eigenvalues of (A,B) are needed, the driver routine\n*  SGEGV should be used instead.  See SGEGV for a description of the\n*  eigenvalues of the generalized nonsymmetric eigenvalue problem\n*  (GNEP).\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVSL  (input) CHARACTER*1\n*          = 'N':  do not compute the left Schur vectors;\n*          = 'V':  compute the left Schur vectors (returned in VSL).\n*\n*  JOBVSR  (input) CHARACTER*1\n*          = 'N':  do not compute the right Schur vectors;\n*          = 'V':  compute the right Schur vectors (returned in VSR).\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VSL, and VSR.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA, N)\n*          On entry, the matrix A.\n*          On exit, the upper quasi-triangular matrix S from the\n*          generalized real Schur factorization.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) REAL array, dimension (LDB, N)\n*          On entry, the matrix B.\n*          On exit, the upper triangular matrix T from the generalized\n*          real Schur factorization.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  ALPHAR  (output) REAL array, dimension (N)\n*          The real parts of each scalar alpha defining an eigenvalue\n*          of GNEP.\n*\n*  ALPHAI  (output) REAL array, dimension (N)\n*          The imaginary parts of each scalar alpha defining an\n*          eigenvalue of GNEP.  If ALPHAI(j) is zero, then the j-th\n*          eigenvalue is real; if positive, then the j-th and (j+1)-st\n*          eigenvalues are a complex conjugate pair, with\n*          ALPHAI(j+1) = -ALPHAI(j).\n*\n*  BETA    (output) REAL array, dimension (N)\n*          The scalars beta that define the eigenvalues of GNEP.\n*          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and\n*          beta = BETA(j) represent the j-th eigenvalue of the matrix\n*          pair (A,B), in one of the forms lambda = alpha/beta or\n*          mu = beta/alpha.  Since either lambda or mu may overflow,\n*          they should not, in general, be computed.\n*\n*  VSL     (output) REAL array, dimension (LDVSL,N)\n*          If JOBVSL = 'V', the matrix of left Schur vectors Q.\n*          Not referenced if JOBVSL = 'N'.\n*\n*  LDVSL   (input) INTEGER\n*          The leading dimension of the matrix VSL. LDVSL >=1, and\n*          if JOBVSL = 'V', LDVSL >= N.\n*\n*  VSR     (output) REAL array, dimension (LDVSR,N)\n*          If JOBVSR = 'V', the matrix of right Schur vectors Z.\n*          Not referenced if JOBVSR = 'N'.\n*\n*  LDVSR   (input) INTEGER\n*          The leading dimension of the matrix VSR. LDVSR >= 1, and\n*          if JOBVSR = 'V', LDVSR >= N.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,4*N).\n*          For good performance, LWORK must generally be larger.\n*          To compute the optimal value of LWORK, call ILAENV to get\n*          blocksizes (for SGEQRF, SORMQR, and SORGQR.)  Then compute:\n*          NB  -- MAX of the blocksizes for SGEQRF, SORMQR, and SORGQR\n*          The optimal LWORK is  2*N + N*(NB+1).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1,...,N:\n*                The QZ iteration failed.  (A,B) are not in Schur\n*                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should\n*                be correct for j=INFO+1,...,N.\n*          > N:  errors that usually indicate LAPACK problems:\n*                =N+1: error return from SGGBAL\n*                =N+2: error return from SGEQRF\n*                =N+3: error return from SORMQR\n*                =N+4: error return from SORGQR\n*                =N+5: error return from SGGHRD\n*                =N+6: error return from SHGEQZ (other than failed\n*                                                iteration)\n*                =N+7: error return from SGGBAK (computing VSL)\n*                =N+8: error return from SGGBAK (computing VSR)\n*                =N+9: error return from SLASCL (various places)\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_jobvsl = argv[0];
  rb_jobvsr = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_lwork = argv[4];

  jobvsl = StringValueCStr(rb_jobvsl)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  lwork = NUM2INT(rb_lwork);
  jobvsr = StringValueCStr(rb_jobvsr)[0];
  ldvsr = lsame_(&jobvsr,"V") ? n : 1;
  ldvsl = lsame_(&jobvsl,"V") ? n : 1;
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
  {
    int shape[2];
    shape[0] = ldvsl;
    shape[1] = n;
    rb_vsl = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vsl = NA_PTR_TYPE(rb_vsl, real*);
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

  sgegs_(&jobvsl, &jobvsr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_alphar, rb_alphai, rb_beta, rb_vsl, rb_vsr, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_sgegs(VALUE mLapack){
  rb_define_module_function(mLapack, "sgegs", rb_sgegs, -1);
}
